Niko Hauzenberger and Michael Pfarrhofer, "Bayesian state-space modeling for analyzing heterogeneous network effects of US monetary policy", 
Scandinavian Journal of Economics, forthcoming.

Michael Weber kindly provided us with the original industry-level dataset and the exogenous monetary policy measure used in 
Ozdagli and Weber (2020), "Monetary policy through production networks: Evidence from the stock market," NBER Working Paper 23424(updated).

The cross-sectional dependency structure via the weighting matrix is based on IO-tables capturing dollar trade flows between industries.
The industries are selected based on the availability of input-output (IO)-tables published
by the Bureau of Economic Analysis (BEA) and the United States Department of Commerce. 
These tables are needed to calculate the cross-sectional linkages in W[t]. 
They are published every five years, and we utilize their 1992, 1997 and 2002 versions.

As exogenous measure of the monetary policy shocks, we rely on high-frequency changes
in Federal funds futures. 
Our information set includes data between February 1994 and December 2008, that is, T = 120. 
The sample starts in 1994 because the Federal Reserve (Fed) changed its communication strategy at this time 
and tick-by-tick stock market data is not available prior to 1993. 
It ends in 2008 to exclude the period when the Fed started its various unconventional monetary policy measures when approaching the zero lower bound.


We provide the readymade dataset as .rda-objects. 

The four objects containing W[t] are: 
- dataTVW.rda: In our baseline model, we allow for time-variation in W[t]. We use the 1997 IO-tables from 1994 to the last FOMC announcement in 2001, 
and rely on the 2002 IO-tables from this point onwards.
- data1992.rda: W[t] = W is based on the 1992 IO-tables.
- data1997.rda: W[t] = W is based on the 1997 IO-tables.
- data2002.rda: W[t] = W is based on the 2002 IO-tables.

The object MPData.rda contains exogenous measure of the monetary policy shocks. 

The R-files contain functions to replicate all results in the paper. Setting run = 1 estimates the proposed most flexible version of the model featuring all types of heterogeneities. Note that the code contains several additional features that were present in earlier versions of the paper.
