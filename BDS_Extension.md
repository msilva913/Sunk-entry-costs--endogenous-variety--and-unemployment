**BFS**

1. Convertion
   
   Monthly rate Gt=（Xt-Xt-1）/Xt-1
   
   year over year monthly growth rate
   
   Index Rate；For example，2006=1

2. Source
   
 Asturias, J., Dinlersoz, E. M., Haltiwanger, J. C., & Hutchinson, R. (2021). Business Applications as Economic Indicators. US Census Bureau, Center for Economic Studies.

 Haltiwanger, J. C. (2022). Entrepreneurship during the COVID-19 pandemic: Evidence from the business formation statistics. Entrepreneurship and Innovation Policy and the Economy, 1(1), 9-42.

 BFS Release[https://www.census.gov/econ/bfs/current/index.html]

 Other...

**BDS**

1. Description
   
   *firms*: A simple count of the number of firms in the cell. For state level tables, a firm with establishmentsin multiple states be counted multiple times, once in each state, irrespective of the portion of the firm residing in that state.
   *emp*:Paid employment consists of full and part-time employees, including salaried officers andexecutives of corporations, who are on the payroll in the pay period including March 12. Includedare employees on paid sick leave, holidays, and vacations.
   *estabs*: A simple count of the number of establishments in the cell.
   *denom*: Davis-Haltiwanger-Schuh (DHS) denominator. For time t, denom is the average of employment for times t and t-1. This variable attempts to prevent transitory shocks from creating a bias to therelationship between net growth from t-1 to t and size.

   *estabs_entry*: A count of establishments born within the cell during the last 12 months.
   *estabs_entry_rate*: 100 * (estabs_entry at time t divided by the average of estabs at t and t-1)
   *estabs_exit*: A count of establishments exiting from within the cell during the last 12 months.
   *estabs_exit_rate*: 00 * (estabs_exit at time t divided by the average of estabs at t and t-1).

   *job_creation*: Count of all employment gains within the cell from expanding and opening establishments between the week of March 12 of the prior year to the current year.
   *job_creation_rate*: 100 * (job_creation / denom)

   *firmdeath_firms*: Count of firms that have exited in their entirety during the period. All establishments owned by the firm must exit to be considered a firm death. This definition of firm death is narrow and strictlyapplied, so that a firm with 100 establishments would not qualify as a firm death if 99 exited while1 continued under different ownership. Note firm legal entities that cease to exist because of merger and acquisition activity are not classified as firm deaths in the BDS data.
   *firmdeath_estabs*: Count of establishments associated with firm deaths.

2. Extension
   
   *firms_exit_rate*: 100 * (firmdeath_firms at time t divided by the average of firms at t and t-1)
   *firms_entry*: firms at time t added with firmdeath_firms at time t minus firms at time t-1. (firms t = firms t-1 + Entry t - Exit t)
   *firms_entry_rate*: 100 * (firms_entry at time t divided by the average of firms at t and t-1)
