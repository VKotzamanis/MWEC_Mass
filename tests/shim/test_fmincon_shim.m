function test_fmincon_shim()
%TEST_FMINCON_SHIM Check the fmincon/optimoptions shims on problems with known solutions.
%   Tolerance 1e-6 on x and multipliers: sqp differentiates by forward differences with step
%   sqrt(eps) = 1.5e-8, so a quadratic objective carries a gradient error of about 1.5e-8 per
%   unit curvature; 1e-6 leaves a margin of more than 50 over that error and is set by the
%   solver, not by the problem.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root, 'tests', 'octave_shims'));
    tol = 1e-6;
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
                        'MaxIterations', 200, 'OptimalityTolerance', 1e-6, ...
                        'ConstraintTolerance', 1e-6, 'StepTolerance', 1e-10);

    % 1. Linear equality + nonlinear inequality + bound, all active.
    %    min (x1-2)^2 + (x2-1)^2 + (x3-3)^2  s.t. x1+x2+x3 = 3, x1^2 <= 1.44, x3 <= 1.
    %    Solution (1.2, 0.8, 1); multipliers (MATLAB sign): eq 0.4, ineqnonlin 0.5, upper(3) 3.6.
    f1 = @(x) (x(1)-2)^2 + (x(2)-1)^2 + (x(3)-3)^2;
    nl1 = @(x) deal(x(1)^2 - 1.44, []);
    [x, fval, ef, out, lam] = fmincon(f1, [0; 0; 0], [], [], [1 1 1], 3, [-5; -5; -5], [5; 5; 1], nl1, opts);
    report('P1 eq+ineq+bound', x, [1.2; 0.8; 1.0], ef, out);
    check('P1 fval', fval, f1([1.2; 0.8; 1]), tol);
    check('P1 lambda.eqlin', lam.eqlin, 0.4, tol);
    check('P1 lambda.ineqnonlin', lam.ineqnonlin, 0.5, tol);
    check('P1 lambda.upper', lam.upper, [0; 0; 3.6], tol);
    check('P1 lambda.lower', lam.lower, [0; 0; 0], tol);
    check('P1 x', x, [1.2; 0.8; 1.0], tol);
    check_converged('P1 exitflag', ef);
    if out.constrviolation > 1e-8
        error('test_fmincon_shim:P1', 'constraint violation %.3e exceeds the 1e-8 request', out.constrviolation);
    end

    % 2. Infeasible start, linear inequality in A*x <= b form.
    %    min (x1-1)^2 + (x2-1)^2 s.t. -x1 - 2 x2 <= -6. Solution (1.6, 2.2), multiplier 1.2.
    f2 = @(x) (x(1)-1)^2 + (x(2)-1)^2;
    [x, ~, ef, out, lam] = fmincon(f2, [-4; -4], [-1 -2], -6, [], [], [], [], [], opts);
    report('P2 infeasible start, A*x<=b', x, [1.6; 2.2], ef, out);
    check('P2 x', x, [1.6; 2.2], tol);
    check('P2 lambda.ineqlin', lam.ineqlin, 1.2, tol);
    check_converged('P2 exitflag', ef);

    % 3. Infeasible start, nonlinear equality.
    %    min x1 + x2 s.t. x1^2 + x2^2 = 2. Solution (-1, -1), multiplier 0.5.
    f3 = @(x) x(1) + x(2);
    nl3 = @(x) deal([], x(1)^2 + x(2)^2 - 2);
    [x, ~, ef, out, lam] = fmincon(f3, [-3; -0.5], [], [], [], [], [], [], nl3, opts);
    report('P3 infeasible start, ceq', x, [-1; -1], ef, out);
    check('P3 x', x, [-1; -1], tol);
    check('P3 lambda.eqnonlin', lam.eqnonlin, 0.5, tol);
    check_converged('P3 exitflag', ef);

    % 4. Bounds only, row-vector x0: output keeps the row shape and fun sees a row.
    %    min (x1-3)^2 + (x2+1)^2 on 0 <= x <= 2. Solution (2, 0), multipliers upper(1) 2, lower(2) 2.
    f4 = @(x) check_row(x) + (x(1)-3)^2 + (x(2)+1)^2;
    [x, ~, ef, out, lam] = fmincon(f4, [1 1], [], [], [], [], [0 0], [2 2], [], opts);
    report('P4 bounds only, row x0', x, [2 0], ef, out);
    if ~isrow(x), error('test_fmincon_shim:shape', 'x must keep the row shape of x0'); end
    check('P4 x', x, [2 0], tol);
    check('P4 lambda.upper', lam.upper, [2; 0], tol);
    check('P4 lambda.lower', lam.lower, [0; 2], tol);
    check_converged('P4 exitflag', ef);

    % 5. Iteration cap: exactly MaxIterations iterations, exitflag 0.
    rosen = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
    capped = optimoptions('fmincon', 'MaxIterations', 3);
    [~, ~, ef, out] = fmincon(rosen, [-1.2; 1], [], [], [], [], [], [], [], capped);
    report('P5 iteration cap', [], [], ef, out);
    check_flag('P5 exitflag', ef, 0);
    check_flag('P5 iterations', out.iterations, 3);

    % 6. OutputFcn gets 'init' (iteration 0) then 'done', and a stop request at 'init' gives -1.
    log = {};
    function_log = @(x, ov, state) log_state(state, ov);
    log_state('reset', []);
    with_fcn = optimoptions('fmincon', 'OutputFcn', function_log);
    [~, ~, ef, out] = fmincon(f2, [-4; -4], [-1 -2], -6, [], [], [], [], [], with_fcn);
    log = log_state('get', []);
    fprintf('P6 OutputFcn states: %s (iterations at done = %d)\n', strjoin(log(:, 1)', ','), log{end, 2});
    if ~isequal(log(:, 1)', {'init', 'done'}) || log{1, 2} ~= 0 || log{end, 2} ~= out.iterations
        error('test_fmincon_shim:outputfcn', 'OutputFcn state sequence or iteration counts wrong');
    end
    stopper = optimoptions('fmincon', 'OutputFcn', @(x, ov, state) true);
    [x, ~, ef] = fmincon(f2, [-4; -4], [-1 -2], -6, [], [], [], [], [], stopper);
    check_flag('P6 stop at init exitflag', ef, -1);
    check('P6 stop at init x', x, [-4; -4], 0);

    % 8. Contradictory constraints x <= 1 and x >= 2: no feasible point, exitflag -2.
    warning('off', 'Octave:SQP-QP-subproblem');
    [~, ~, ef, out] = fmincon(@(x) x^2, 0, [1; -1], [1; -2], [], [], [], [], [], opts);
    warning('on', 'Octave:SQP-QP-subproblem');
    report('P8 infeasible problem', [], [], ef, out);
    check_flag('P8 exitflag', ef, -2);
    if out.constrviolation < 0.5
        error('test_fmincon_shim:P8', 'a violation of at least 0.5 is unavoidable, got %.3e', out.constrviolation);
    end

    % 7. optimoptions keeps the solver name and updates an existing struct.
    o = optimoptions('fmincon', 'MaxIterations', 5);
    o = optimoptions(o, 'MaxIterations', 7, 'Display', 'off');
    check_flag('P7 SolverName', strcmp(o.SolverName, 'fmincon'), true);
    check_flag('P7 updated option', o.MaxIterations, 7);
end

function v = check_row(x)
    if ~isrow(x), error('test_fmincon_shim:shape', 'fun received a non-row x'); end
    v = 0;
end

function varargout = log_state(state, ov)
    persistent entries
    switch state
        case 'reset'
            entries = {};
            return;
        case 'get'
            varargout{1} = entries;
            return;
    end
    entries(end+1, :) = {state, ov.iteration};
    varargout{1} = false;
end

function report(name, x, x_exact, ef, out)
    if isempty(x)
        fprintf('%s: exitflag %d, iterations %d, funcCount %d, constrviolation %.2e\n', ...
                name, ef, out.iterations, out.funcCount, out.constrviolation);
    else
        fprintf('%s: exitflag %d, iterations %d, funcCount %d, constrviolation %.2e, max|x-x*| %.2e\n', ...
                name, ef, out.iterations, out.funcCount, out.constrviolation, max(abs(x(:) - x_exact(:))));
    end
end

function check(name, value, expected, tol)
    err = max(abs(value(:) - expected(:)));
    fprintf('  %-26s max error %.3e\n', name, err);
    if err > tol
        error('test_fmincon_shim:value', '%s: error %.3e exceeds %.1e', name, err, tol);
    end
end

function check_converged(name, ef)
% sqp reports 101 (exitflag 1) or, when the step stalls at the solution, 104 (exitflag 2).
    fprintf('  %-26s %d\n', name, ef);
    if ~any(ef == [1 2])
        error('test_fmincon_shim:flag', '%s: got %d, expected 1 or 2', name, ef);
    end
end

function check_flag(name, value, expected)
    if value ~= expected
        error('test_fmincon_shim:flag', '%s: got %g, expected %g', name, value, expected);
    end
end
