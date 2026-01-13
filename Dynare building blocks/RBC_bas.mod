/*
 * This file presents a baseline RBC model with TFP and government spending shocks, calibrated to US data from
 *  1947Q4:2016Q1. The model setup is described in Handout_RBC_model.pdf and resembles the one in King/Rebelo (1999): 
 *  Resuscitating Real Business Cycles, Handbook of Macroeconomics, Volume 1 and
 *  Romer (2012), Advanced macroeconomics, 4th edition
 *  The driving processes are estimated as AR(1)-processes on linearly detrended data.
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2016 Johannes Pfeifer
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
 */

//****************************************************************************
//Define variables
//****************************************************************************

var y           ${y}$ (long_name='output')
    c           ${c}$ (long_name='consumption')
    k           ${k}$ (long_name='capital')
    l           ${l}$ (long_name='hours')
    z           ${z}$ (long_name='TFP')
    r           ${r}$ (long_name='annualized interest rate')
    w           ${w}$ (long_name='real wage')
    invest      ${i}$ (long_name='investment') 
    log_y       ${\log(y)}$ (long_name='log output')
    log_k       ${\log(k)}$ (long_name='log capital stock')
    log_c       ${\log(c)}$ (long_name='log consumption')
    log_l       ${\log(l)}$ (long_name='log labor')
    log_w       ${\log(w)}$ (long_name='log real wage')
    log_invest  ${\log(i)}$ (long_name='log investment')
    ;

varexo eps_z ${\varepsilon_z}$ (long_name='TFP shock')
    ;
    
parameters 
    beta    ${\beta}$   (long_name='discount factor')
    chi     ${\chi}$    (long_name='labor disutility parameter')
    psi     ${\psi}$    (long_name='Frisch elasticity of labor supply')
    sigma   ${\sigma}$  (long_name='risk aversion')
    delta   ${\delta}$  (long_name='depreciation rate')
    alpha   ${\alpha}$  (long_name='capital share')
    rhoz    ${\rho_z}$  (long_name='persistence TFP shock')
    ;
//****************************************************************************
//Set parameter values
//****************************************************************************

beta = 0.99;
delta = 0.025;
chi = 0.9241;
psi = 4.0;
sigma=1;                // risk aversion
alpha= 0.33;            // capital share
rhoz=0.979;              //technology autocorrelation base on linearly detrended Solow residual

//****************************************************************************
//enter the model equations (model-block)
//****************************************************************************

model;
[name='Euler equation']
c^(-sigma)=beta*c(+1)^(-sigma)*(alpha*exp(z(+1))*(k/l(+1))^(alpha-1)+(1-delta));

[name='Labor FOC']
chi*l^(1/psi)=w*c^(-sigma);

[name='Law of motion capital'] 
k=(1-delta)*k(-1)+invest;

[name='resource constraint']
y=invest+c;

[name='production function']
y=exp(z)*k(-1)^alpha*l^(1-alpha);

[name='real wage/firm FOC labor']
w=(1-alpha)*y/l;

[name='annualized real interest rate/firm FOC capital']
r=4*alpha*y/k(-1);

[name='exogenous TFP process']
z=rhoz*z(-1)+eps_z;

[name='Definition log output']
log_y = 100*log(y);

[name='Definition log capital']
log_k = 100*log(k);

[name='Definition log consumption']
log_c = 100*log(c);

[name='Definition log hours']
log_l = 100*log(l);

[name='Definition log wage']
log_w = 100*log(w);

[name='Definition log investment']
log_invest = 100*log(invest);
end;

//****************************************************************************
// Provide steady state values and calibrate the model to steady state labor of 0.33, 
// i.e. compute the corresponding steady state values
// and the labor disutility parameter by hand;
//****************************************************************************

initval;
    //Do Calibration
    l=1.0;
    k = ((1/beta-(1-delta))/alpha)^(1/(alpha-1))*l; 
    invest = delta*k;

    y=k^alpha*l^(1-alpha);
    c = y-invest;
    w = (1-alpha)*y/l;
    r = 4*alpha*y/k;
    log_y = 100*log(y);
    log_k = 100*log(k);
    log_c = 100*log(c);
    log_l = 100*log(l);
    log_w = 100*log(w);
    log_invest = log(invest);
    z = 0; 
end;

//****************************************************************************
//set shock variances
//****************************************************************************

shocks;
    //var eps_z=0.66^2;
    var eps_z = 0.01^2;
end;

//****************************************************************************
//check the starting values for the steady state
//****************************************************************************

resid;

//****************************************************************************
// compute steady state given the starting values
//****************************************************************************

steady;
//****************************************************************************
// check Blanchard-Kahn-conditions
//****************************************************************************

check;

//****************************************************************************
// compute policy function at first order, do IRFs and compute moments with HP-filter
//****************************************************************************

stoch_simul(order=1,irf=40) log_y  log_l  z;