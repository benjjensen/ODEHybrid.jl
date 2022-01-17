function a = test_ode(t, x, u, tsl)
    tsl.add('A', t, x, true, 'group1');
    a = ([0 1; 2 0] * x) + [0; 1] * u{:};
end