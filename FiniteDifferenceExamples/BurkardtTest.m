fd1d_wave_test();

function fd1d_wave_test02 ( )

%*****************************************************************************80
%
%% FD1D_WAVE_TEST02 tests the FD1D finite difference wave computation.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 January 2012
%
%  Author:
%
%    John Burkardt
%
  x_num = 16;
  x1 = 0.0;
  x2 = 1.5;
  x_vec = linspace ( x1, x2, x_num );

  t_num = 41;
  t1 = 0.0;
  t2 = 4.0;
%
%  Changing T2 to 4.5 is enough to push the algorithm into instability.
%
%  t2 = 4.5;
%
  t_vec = linspace ( t1, t2, t_num );
  t_delta = ( t2 - t1 ) / ( t_num - 1 );

  c = 1.0;
  alpha = fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c );

  u = zeros ( t_num, x_num );
%
%  Load the initial condition.
%
  u1(1:x_num) = u_t1_02 ( x_num, x_vec );
  u(1,1:x_num) = u1(1:x_num);
%
%  Take the first step.
%
  t = t_vec(2);
  u2(1:x_num) = fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, @u_x1_02, ...
    @u_x2_02, @ut_t1_02, u1 );
  u(2,1:x_num) = u2(1:x_num);
%
%  Take all the other steps.
%
  for i = 3 : t_num
    t = t_vec(i);
    u3(1:x_num) = fd1d_wave_step ( x_num, t, alpha, @u_x1_02, @u_x2_02, u1, u2 );
    u(i,1:x_num) = u3(1:x_num);
    u1(1:x_num) = u2(1:x_num);
    u2(1:x_num) = u3(1:x_num);
  end

  u_min = min ( min ( u ) );
  u_max = max ( max ( u ) );
%
%  Plot the solution as it evolves in time.
%
  x = linspace ( x1, x2, x_num );
  for t = 1 : t_num
    time = ( ( t_num - t ) * t1 + ( t - 1 ) * t2 ) / ( t_num - 1 );
    plot ( x(1:x_num), u(t,1:x_num), 'b*- ' )
    grid on
    axis ( [ x1, x2, u_min-1.0, u_max+1.0 ] );
    title ( sprintf ( 'Step %d, Time %f\n', t - 1, time ) )
    xlabel ( '<-- X -->' )
    ylabel ( 'Vertical displacement' )
    pause
  end
%
%  Plot the entire solution as a surface.
%
  if ( 1 )

    t_vec = linspace ( t1, t2, t_num );
    x_vec = linspace ( x1, x2, x_num );
    
    [ T, X ] = meshgrid ( t_vec, x_vec );

    surfl ( T', X', u )

    axis ( [ t1, t2, x1, x2, u_min, u_max ] );
    xlabel ( '<-- T -->' );
    ylabel ( '<-- X -->' );
    zlabel ( 'Vertical displacement' );
    title ( 'Solution U(T,X)' );

  end
%
%  Write the solution to a file.
%
  r8mat_write ( 'test02_plot.txt', t_num, x_num, u );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Plot data written to "test02_plot.txt".\n' );

  return
end
function u = u_x1_02 ( t )

%*****************************************************************************80
%
%% U_X1_02 evaluates U at the boundary X1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real T, the time.
%
%    Output, real U, the value of U(T,X1).
%
  u = 0.0;

  return
end
function u = u_x2_02 ( t )

%*****************************************************************************80
%
%% U_X2_02 evaluates U at the boundary X2.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real T, the time.
%
%    Output, real U, the value of U(T,X2).
%
  u = 0.0;

  return
end
function u = u_t1_02 ( x_num, x_vec )

%*****************************************************************************80
%
%% U_T1_02 evaluates U at the initial time T1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer X_NUM, the number of nodes.
%
%    Input, real X_VEC(X_NUM), the spatial node coordinates.
%
%    Output, real U(1,X_NUM), the value of U at the initial time,
%    and every node.
%
  u = zeros(1,x_num);
  u(1,1:x_num) = sin ( 2 * pi * x_vec );

  return
end
function ut = ut_t1_02 ( x_num, x_vec )

%*****************************************************************************80
%
%% UT_T1_02 evaluates dUdT at the initial time T1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer X_NUM, the number of nodes.
%
%    Input, real X_VEC(X_NUM), the spatial node coordinates.
%
%    Output, real UT(1,X_NUM), the value of dU/dt at the initial time,
%    and every node.
%
  ut = zeros(1,x_num);

  return
end
function fd1d_wave_test01 ( )

%*****************************************************************************80
%
%% FD1D_WAVE_TEST01 tests the FD1D finite difference wave computation.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
  x_num = 16;
  x1 = 0.0;
  x2 = 1.5;
  x_vec = linspace ( x1, x2, x_num );

  t_num = 41;
  t1 = 0.0;
  t2 = 4.0;
  t_vec = linspace ( t1, t2, t_num );
  t_delta = ( t2 - t1 ) / ( t_num - 1 );

  c = 1.0;
  alpha = fd1d_wave_alpha( x_num, x1, x2, t_num, t1, t2, c );

  u = zeros ( t_num, x_num );
%
%  Load the initial condition.
%
  u1(1:x_num) = u_t1_01 ( x_num, x_vec );
  u(1,1:x_num) = u1(1:x_num);
%
%  Take the first step.
%
  t = t_vec(2);
  u2(1:x_num) = fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, @u_x1_01, ...
    @u_x2_01, @ut_t1_01, u1 );
  u(2,1:x_num) = u2(1:x_num);
%
%  Take all the other steps.
%
  for i = 3 : t_num
    t = t_vec(i);
    u3(1:x_num) = fd1d_wave_step ( x_num, t, alpha, @u_x1_01, @u_x2_01, u1, u2 );
    u(i,1:x_num) = u3(1:x_num);
    u1(1:x_num) = u2(1:x_num);
    u2(1:x_num) = u3(1:x_num);
  end

  u_min = min ( min ( u ) );
  u_max = max ( max ( u ) );
%
%  Plot the solution as it evolves in time.
%
  x = linspace ( x1, x2, x_num );
  for t = 1 : t_num
    time = ( ( t_num - t ) * t1 + ( t - 1 ) * t2 ) / ( t_num - 1 );
    plot ( x(1:x_num), u(t,1:x_num), 'b*- ' )
    grid on
    axis ( [ x1, x2, u_min-1.0, u_max+1.0 ] );
    title ( sprintf ( 'Step %d, Time %f\n', t - 1, time ) )
    xlabel ( '<-- X -->' )
    ylabel ( 'Vertical displacement' )
    pause
  end
%
%  Plot the entire solution as a surface.
%
  if ( 1 )

    t_vec = linspace ( t1, t2, t_num );
    x_vec = linspace ( x1, x2, x_num );
    
    [ T, X ] = meshgrid ( t_vec, x_vec );

    surfl ( T', X', u )

    axis ( [ t1, t2, x1, x2, u_min, u_max ] );
    xlabel ( '<-- T -->' );
    ylabel ( '<-- X -->' );
    zlabel ( 'Vertical displacement' );
    title ( 'Solution U(T,X)' );

  end
%
%  Write the solution to a file.
%
  r8mat_write ( 'test01_plot.txt', t_num, x_num, u );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Plot data written to "test01_plot.txt".\n' );

  return
end
function u = u_x1_01 ( t )

%*****************************************************************************80
%
%% U_X1_01 evaluates U at the boundary X1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real T, the time.
%
%    Output, real U, the value of U(T,X1).
%
  if ( t == 0.10 )
    u = 2.0;
  elseif ( t == 0.20 ) 
    u = 10.0;
  elseif ( t == 0.30 )
    u = 8.0;
  elseif ( t == 0.40 )
    u = 5.0;
  else
    u = 0.0;
  end

  return
end
function u = u_x2_01 ( t )

%*****************************************************************************80
%
%% U_X2_01 evaluates U at the boundary X2.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real T, the time.
%
%    Output, real U, the value of U(T,X2).
%
  u = 0.0;

  return
end
function u = u_t1_01 ( x_num, x_vec )

%*****************************************************************************80
%
%% U_T1_01 evaluates U at the initial time T1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer X_NUM, the number of nodes.
%
%    Input, real X_VEC(X_NUM), the coordinates of the nodes.
%
%    Output, real U(1,X_NUM), the value of U at the initial time,
%    and every node.
%
  u = zeros(1,x_num);

  return
end
function ut = ut_t1_01 ( x_num, x_vec )

%*****************************************************************************80
%
%% UT_T1_01 evaluates dUdT at the initial time T1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer X_NUM, the number of nodes.
%
%    Input, real X_VEC(X_NUM), the coordinates of the nodes.
%
%    Output, real UT(1,X_NUM), the value of dU/dt at the initial time,
%    and every node.
%
  ut = zeros(1,x_num);

  return
end
function fd1d_wave_test ( )

%*****************************************************************************80
%
%% FD1D_WAVE_TEST tests FD1D_WAVE.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    15 January 2019
%
%  Author:
%
%    John Burkardt
%
  addpath ( '../fd1d_wave' );

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'FD1D_WAVE_TEST\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Test FD1D_WAVE.\n' );

  fd1d_wave_test01 ( );
  fd1d_wave_test02 ( );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'FD1D_WAVE_TEST\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../fd1d_wave' );

  return
end
function alpha = fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c )

%*****************************************************************************80
%
%% FD1D_WAVE_ALPHA computes ALPHA for the 1D wave equation.
%
%  Discussion:
%
%    The explicit timestepping procedure uses the quantity ALPHA, which
%    is determined by this function.
%
%    If the spatial region bounds are X1 <= X <= X2, containing X_NUM equally
%    spaced nodes, including the endpoints, and the time domain similarly
%    extends from T1 <= T <= T2 containing T_NUM equally spaced time values,
%    then
%
%      ALPHA = C * DT / DX
%            = C * ( ( T2 - T1 ) / ( T_NUM - 1 ) )
%                / ( ( X2 - X1 ) / ( X_NUM - 1 ) ).
%
%    For a stable computation, it must be the case that ALPHA < 1.
%
%    If ALPHA is greater than 1, then the middle coefficient 1-C^2 DT^2 / DX^2 
%    is negative, and the sum of the magnitudes of the three coefficients becomes 
%    unbounded.  In such a case, the user must reduce the time step size 
%    appropriately.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    George Lindfield, John Penny,
%    Numerical Methods Using MATLAB,
%    Second Edition,
%    Prentice Hall, 1999,
%    ISBN: 0-13-012641-1,
%    LC: QA297.P45.
%
%  Parameters:
%
%    Input, integer X_NUM, the number of nodes in the X direction.
%
%    Input, real X1, X2, the first and last X coordinates.
%
%    Input, integer T_NUM, the number of time steps, including the 
%    initial condition.
%
%    Input, real T1, T2, the first and last T coordinates.
%
%    Input, real C, a parameter which gives the speed of waves.
%
%    Output, real ALPHA, the stability coefficient.
%
  t_delta = ( t2 - t1 ) / ( t_num - 1 );
  x_delta = ( x2 - x1 ) / ( x_num - 1 );
  alpha = c * t_delta / x_delta;

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Stability condition ALPHA = C * DT / DX = %f\n', alpha );

  if ( 1.0 < abs ( alpha ) )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'FD1D_WAVE_ALPHA - Warning!\n' );
    fprintf ( 1, '  The stability condition |ALPHA| <= 1 fails.\n' );
    fprintf ( 1, '  Computed results are liable to be inaccurate.\n' );
  end

  return
end
function u2 = fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, u_x1, u_x2, ...
  ut_t1, u1 )

%*****************************************************************************80
%
%% FD1D_WAVE_START takes the first step for the wave equation.
%
%  Discussion:
%
%    This program solves the 1D wave equation of the form:
%
%      Utt = c^2 Uxx
%
%    over the spatial interval [X1,X2] and time interval [T1,T2],
%    with initial conditions:
%
%      U(T1,X)  = U_T1(X),
%      Ut(T1,X) = UT_T1(X),
%
%    and boundary conditions of Dirichlet type:
%
%      U(T,X1) = U_X1(T),
%      U(T,X2) = U_X2(T).
%
%    The value C represents the propagation speed of waves.
%
%    The program uses the finite difference method, and marches
%    forward in time, solving for all the values of U at the next
%    time step by using the values known at the previous two time steps.
%
%    Central differences may be used to approximate both the time
%    and space derivatives in the original differential equation.
%
%    Thus, assuming we have available the approximated values of U
%    at the current and previous times, we may write a discretized
%    version of the wave equation as follows:
%
%      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
%      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
%
%    If we multiply the first term by C^2 and solve for the single
%    unknown value U(T+dt,X), we have:
%
%      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
%                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
%                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
%                  -                                  U(T-dT,X   )
%
%    (Equation to advance from time T to time T+dT, except for FIRST step!)
%
%    However, on the very first step, we only have the values of U
%    for the initial time, but not for the previous time step.
%    In that case, we use the initial condition information for dUdT
%    which can be approximated by a central difference that involves
%    U(T+dT,X) and U(T-dT,X):
%
%      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
%
%    and so we can estimate U(T-dT,X) as
%
%      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
%
%    If we replace the "missing" value of U(T-dT,X) by the known values
%    on the right hand side, we now have U(T+dT,X) on both sides of the
%    equation, so we have to rearrange to get the formula we use
%    for just the first time step:
%
%      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
%                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
%                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
%                  +  dT *                         dU/dT(T,   X   )
%
%    (Equation to advance from time T to time T+dT for FIRST step.)
%
%    It should be clear now that the quantity ALPHA = C * DT / DX will affect
%    the stability of the calculation.  If it is greater than 1, then
%    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
%    sum of the magnitudes of the three coefficients becomes unbounded.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    George Lindfield, John Penny,
%    Numerical Methods Using MATLAB,
%    Second Edition,
%    Prentice Hall, 1999,
%    ISBN: 0-13-012641-1,
%    LC: QA297.P45.
%
%  Parameters:
%
%    Input, integer X_NUM, the number of nodes in the X direction.
%
%    Input, real X_VEC(X_NUM), the spatial coordinates of the nodes.
%
%    Input, real T, the time after the first step has been taken.
%    In other words, T = T1 + T_DELTA.
%
%    Input, real T_DELTA, the time step.
%
%    Input, real ALPHA, the stability coefficient, computed by FD1D_WAVE_ALPHA.
%
%    Input, real U_X1(T), U_X2(T), functions for the left and right boundary 
%    conditions.
%
%    Input, real UT_T1(X), the function that evaluates dUdT at the initial time.
%
%    Input, real U1(X_NUM), the initial condition.
%
%    Output, real U2(X_NUM), the solution at the first time step.
%
  ut = ut_t1 ( x_num, x_vec );

  u2 = zeros ( 1, x_num );

  u2(1) = u_x1 ( t );
  u2(2:x_num-1) =         alpha^2   * u1(3:x_num) / 2 ...
                  + ( 1 - alpha^2 ) * u1(2:x_num-1) ...
                  +       alpha^2   * u1(1:x_num-2) / 2 ...
                  +        t_delta  * ut(2:x_num-1);
  u2(x_num) = u_x2 ( t );

  return
end
function u3 = fd1d_wave_step ( x_num, t, alpha, u_x1, u_x2, u1, u2 )

%*****************************************************************************80
%
%% FD1D_WAVE_STEP computes a step of the 1D wave equation.
%
%  Discussion:
%
%    This program solves the 1D wave equation of the form:
%
%      Utt = c^2 Uxx
%
%    over the spatial interval [X1,X2] and time interval [T1,T2],
%    with initial conditions:
%
%      U(T1,X)  = U_T1(X),
%      Ut(T1,X) = UT_T1(X),
%
%    and boundary conditions of Dirichlet type:
%
%      U(T,X1) = U_X1(T),
%      U(T,X2) = U_X2(T).
%
%    The value C represents the propagation speed of waves.
%
%    The program uses the finite difference method, and marches
%    forward in time, solving for all the values of U at the next
%    time step by using the values known at the previous two time steps.
%
%    Central differences may be used to approximate both the time
%    and space derivatives in the original differential equation.
%
%    Thus, assuming we have available the approximated values of U
%    at the current and previous times, we may write a discretized
%    version of the wave equation as follows:
%
%      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
%      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
%
%    If we multiply the first term by C^2 and solve for the single
%    unknown value U(T+dt,X), we have:
%
%      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
%                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
%                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
%                  -                                  U(T-dT,X   )
%
%    (Equation to advance from time T to time T+dT, except for FIRST step!)
%
%    However, on the very first step, we only have the values of U
%    for the initial time, but not for the previous time step.
%    In that case, we use the initial condition information for dUdT
%    which can be approximated by a central difference that involves
%    U(T+dT,X) and U(T-dT,X):
%
%      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
%
%    and so we can estimate U(T-dT,X) as
%
%      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
%
%    If we replace the "missing" value of U(T-dT,X) by the known values
%    on the right hand side, we now have U(T+dT,X) on both sides of the
%    equation, so we have to rearrange to get the formula we use
%    for just the first time step:
%
%      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
%                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
%                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
%                  +  dT *                         dU/dT(T,   X   )
%
%    (Equation to advance from time T to time T+dT for FIRST step.)
%
%    It should be clear now that the quantity ALPHA = C * DT / DX will affect
%    the stability of the calculation.  If it is greater than 1, then
%    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
%    sum of the magnitudes of the three coefficients becomes unbounded.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    George Lindfield, John Penny,
%    Numerical Methods Using MATLAB,
%    Second Edition,
%    Prentice Hall, 1999,
%    ISBN: 0-13-012641-1,
%    LC: QA297.P45.
%
%  Parameters:
%
%    Input, integer X_NUM, the number of nodes in the X direction.
%
%    Input, real T, the new time, that is, the current time + T_DELTA.
%
%    Input, real ALPHA, the stability coefficient, computed by FD1D_WAVE_ALPHA.
%
%    Input, real U_X1(T), U_X2(T), functions for the left and right boundary 
%    conditions.
%
%    Input, real U1(X_NUM), the solution at the old time.
%
%    Input, real U2(X_NUM), the solution at the current time.
%
%    Output, real U3(X_NUM), the solution at the new time.
%
  u3 = zeros ( 1, x_num );

  u3(1)         = u_x1 ( t );
  u3(2:x_num-1) =             alpha^2   * u2(3:x_num) ...
                  + 2 * ( 1 - alpha^2 ) * u2(2:x_num-1) ...
                  +           alpha^2   * u2(1:x_num-2) ...
                  -                       u1(2:x_num-1);
  u3(x_num)     = u_x2 ( t );

  return
end
function yv = piecewise_linear ( nd, xd, yd, nv, xv )

%*****************************************************************************80
%
%% PIECEWISE_LINEAR evaluates a piecewise linear spline.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    24 January 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer ND, the number of data points.
%
%    Input, real XD(ND), YD(ND), the data values.
%
%    Input, integer NV, the number of evaluation points.
%
%    Input, real XV(NV), the evaluation arguments.
%
%    Output, real YV(NV), the values.
%
  yv = zeros ( nv );

  for iv = 1 : nv

    if ( xv(iv) < xd(1) )
      yv(iv) = yd(1);
    elseif ( xd(nd) < xv(iv) )
      yv(iv) = yd(nd);
    else

      for id = 2 : nd
        if ( xv(iv) < xd(id) )
          yv(iv) = ( ( xd(id) - xv(iv)            ) * yd(id-1) ...
                   + (          xv(iv) - xd(id-1) ) * yd(id) ) ...
                   / ( xd(id)          - xd(id-1) );
          break
        end
      end

    end

  end

  return
end
function r8mat_write ( output_filename, m, n, table )

%*****************************************************************************80
%
%% R8MAT_WRITE writes an R8MAT file.
%
%  Discussion:
%
%    An R8MAT is an array of R8's.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 February 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, string OUTPUT_FILENAME, the output filename.
%
%    Input, integer M, the spatial dimension.
%
%    Input, integer N, the number of points.
%
%    Input, real TABLE(M,N), the points.
%

%
%  Open the file.
%
  output_unit = fopen ( output_filename, 'wt' );

  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8MAT_WRITE - Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'R8MAT_WRITE - Error!' );
  end
%
%  Write the data.
%
%  For smaller data files, and less precision, try:
%
%     fprintf ( output_unit, '  %14.6e', table(i,j) );
%
  for j = 1 : n
    for i = 1 : m
      fprintf ( output_unit, '  %24.16e', table(i,j) );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  Close the file.
%
  fclose ( output_unit );

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end
