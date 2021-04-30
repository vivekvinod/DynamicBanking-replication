# DynamicBanking-replication

Replicating the numerical simulations of Dynamic Banking and the Value of Deposits (2020)  

Dynamicbanking.jl shows the final simplified code that attempts to solve the model.

Bolton2020.ipynb - Plotting the graphs in the report

https://github.com/matthieugomez/EconPDEs.jl ECONPDEs for the package and sample code that we edited.

We get the sharp transitions as plotted in the original paper when we optimize the Boundary condition 4 (out[4]) by giving the values of the function at the dividend payout boundary. In the updated version of the code, we can find that the marginal value of Equity shoots up to a value of 12 ( slightly greater than Bolton (2020) )at the equity issuance boundary.
