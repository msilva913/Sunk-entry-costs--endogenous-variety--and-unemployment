# Program files for Sunk vacancy costs, endogenous product variety, and labor market frictions 
## Key contribution
We show how to
- construct the empirical series and obtain summaries and plots
- calculate the steady state, calibrate the model, perform comparative statics, and generate tables
- solve the model via perturbation and generate impulse responses and simulations
- compare impulse responses with appropriate recalibration
- **To do: estimate by simulated method of moments**

## File documentation.
We document the files into two types, those necessary for descriptive statistics and those necessary for model solution and estimation
### Descriptive statistics
- Main file is `observables.py`. 
	- The function `construct_data` constructs the time series for consumption, consumptions share, unemployment, vacancies, job finding rate, labor productivity, labor share, separation rate, and wages. Includes
		- Plot of unemployment and vacancies
		- Beveridge curve
		- Estimated elasticity of Beveridge curve
		- Estimated matching function
		- Plot of unemployment and labor productivity
	- The function `create_stats_table(data, caption, label)` generates table of second moments
	- Function `first_momens(dat, alpha_hat)` summarizes data means.
	- Calls `time_series_functions.py` for calculation of moments and filtering
	- `BFS.py` provides additional plots and descriptive stats of Business Formation Statistics
### Analysis
#### Steady state and calibration
- `steady_state.jl` provides main functions characterizing the steady state, the calibration file, and also writes a calibration table
	- Function `θ_fun` solves for steady-state theta and serves as core part of `steady_state`
	- Function `calibrate_labor_share` does main calibration, including labor share of income as target
	- There is also `calibrate` function, which does not include labor share  of income as target
	- Functions `N_jcc(θ, para)` and `N_res(θ, para)` are the ss job creation and resource constraint curves
	- `calibration_table` prepares calibration table for paper and slides
- `steady_state_checks.jl` does various accuracy checks regarding steady state and calibration
- `steady_state_flow.jl` defines the steady state of model with flow vacancy posting costs as in DMP
- `generate_tables.jl` creates table of steady state shares given the calibrated steady state
#### Solution and simulation
- `run_solution.jl` is the main program to solve the model via perturbation
	- Builds on the interface by Alvaro Salazar-Perez and Hernán D. Seoane
	- `solution_functions.jl` contains library of functions to process and solve the model
	- Idea
		- Use symbolic variables to declare parameters, dependent parameters and variables
		- Parameters include structural parameters and shock processes
		- Write variables explicitly as states $x$ and jumpers $y$, together with future values $xp$ and $yp$ 
		- Write array $f$ of model equations in form $E(x,y,xp,yp) = 0$ 
		- Write array of steady-state values 
		- Put information togther in `model` named tuple 
		- Run model through processor. This interprets the symbolic variables and also automatically expresses things in log deviations if specified
		- Use `solution_interface` to solve model (key function is `solve_model`)
	- Use `simulate_model` to both simulate artificial data and generate impulse responses.
	- `Impulse_response_plots.jl` and `Impulse_response_comparison.jl` provide main plotting functions
		- For IRF's, choose `flag_IR=true` and `flag_logdev=true` 
	- Save output using `serialize`
- `Second_moments.jl' using saved output to calculate moments (ie, HP-filtered)
- Makes use of functions `from time_series_fun.jl`
