DeepSeaVents README

This project was inspired by impermanent habitat patches that support metapopulations of species (e.g. deep sea hydrothermal vents) and provides code for calculating metrics of metapopulation persistence and size and simulating metapopulation dynamics as the configuration of habitable patches changes through individual patches randomly and independently turning on and off at any given rates. 

## Things you might want to do and recipes for how to use this code to do them:
### A) Make a landscape specifying a spatial structure of permanent or impermanent habitat patches
Simply run the make_landscape( ) function with parameters set as desired for how you would like the landscape to be spatially structured (see #make_landscape function in the functions section for details).
### B) Calculate metrics of metapopulation size and persistence in a landscape of permanent habitat
If you have not already done so make a **landscape** describing the spatial structure of habitat patches in the landscape (see recipe A). Then you can just run the **get_MetapopMetrics( )** with parameters set as desired (see #get_MetapopMetrics in the functions section for details). No need to do any of the following steps, they are done for you within this function.
#### B1) Calculate a matrix of inter patch distances in a landscape
This is actually already done for you by the make_landscape( ) function but you can also use the **get_distmat( )** function to slice off just the distance matrix from the **landscape** dataframe made by the make_landscape( ) function if you just want this for something.
#### B2) Calculate the dispersal kernel for a species in landscape
You just need a distance matrix of interpatch distances (**dist.mat**) to supply to the **get_dispkernel( )** (see recipe B1 if you haven’t got this) and can then set the other parameters as desired for how the species should disperse in a landscape (see the #get_dispkernel in the functions section for details).  
#### B3) Calculate the landscape capacity (*lambda.M*)
You just need a **landscape** as created by the **make_landscape( )** function (see A) and a dispersal kernel (**disp.kernel**) as calculated by the **get_dispkernel( )** function (see B2) for that species in that landscape (see the make_landscape and get_dispkernel functions in the functions section for details). Supply these to the **get_lambdaM( )** function along with the other species specific parameters of this function set as desired (see #get_lambdaM in the functions section for details)
#### B4) Calculate persistence capacity (*lambda.M.delta*)
This is just the landscape capacity (**lambda.M**) x **delta** (**e.rate**/**c.rate**) so all you need to do is follow recipe B3 and multiply lambda.M accordingly.
#### 5) Calculate equilibrium occupancy (*pstar*) of a metapopulation
You just need a **landscape** as created by the **make_landscape( )** function (see A) and a dispersal kernel (**disp.kernel**) as calculated by the **get_dispkernel( )** function (see B2) for that species in that landscape (see the make_landscape and get_dispkernel functions in the functions section for details). Supply these to the **get_pstar( )** function along with the other species specific parameters of this function and the parameter for how many times to iterate the function for finding pstar (more iterations more accurate) set as desired (see #get_pstar in the functions section for details)
### C) Solve the SRLM for the metapopulation dynamics over a period of time in a static landscape of permanent patches
If you have not already done so make a **landscape** describing the spatial structure of habitat patches in the landscape (see recipe A), you can just run the #setNrun_SRLMODE with parameters set as desired (see #setNrun_SRLMODE in the functions section for details). 
### D) Create the generator matrix of a CTMC describing the instantaneous transition probabilities between spatial configurations of habitat when patches are impermanent
First consider the number of impermanent patches that may be in the system. If this number is small (<15-20) feel free to use the **make_Gmat( )** to describe the instantaneous transition probabilities from and to every single spatial configuration the habitat could have, otherwise you will need to use the **make_reducedGmat( )** and will only be able to describe the instantaneous transition probabilities between how many patches may be present in habitat configurations. Don’t worry, the latter can also be used to obtain many properties of the former.
### E) Calculate the QED of spatial configurations of habitat when patches are impermanent
You only need the submatrix of the generator matrix describing the instantaneous probabilities to and from the non-absorbing states (aka just excluding the no habitat state). Just supply the generator matrix, **G.mat** , as made by recipe D to the **Gmat_to_Csubmat( )** function and supply the output (**C.submat**) to the **get_QED( )** function. Importantly, the output **QED** will either describe the probabilities of each spatial configuration (if you used the **make_Gmat( )** function) or of just how many patches are in a spatial configuration (if you used the **make_reduced_Gmat( )** function) of the landscape at QED. The probability of each spatial configuration with the same number of patches present in it is the same at QED since patch impermanence is assumed to be random and independent. Therefore, you only need to divide the probabilities of however many patches may be present in configurations at QED by the number of configurations possible with only that many patches present to convert these to probabilities of each particular spatial configuration at QED identified by how many patches are present in it if you used the **make_reduced_Gmat**. The **reduced_to_full_QED( )** function will do that or you (see #reduced_to_full_QED in the functions section).
### F) Calculate the average time till no habitat is present when patches are impermanent
You only need the submatrix of the generator matrix describing the instaneous probabilities to and from the non-absorbing states (aka just excluding the no habitat state). Just supply the generator matrix, **G.mat** , as made by recipe D to the **Gmat_to_Csubmat( )** function and supply the output (**C.submat**) to the **get_Tabsorp( )** function. The output vector gives the times to absorption (the no-habitat state) from each of the other states the generator matrix describes (whether these be representative of specific spatial configurations of patches or simply how many patches present within them) (refer to whether you used **make_Gmat( )** or **make_reduced_Gmat( )** in recipe D).
### G) Check the rate of convergence to QED will be fast relative to the rate of convergence to the absorbing (no-habitat) state
Follow recipe E, the second item in the output list provides a metric (**rho1vrho2x2**) suggesting this is the case when it is <1 (see the Glossary of Function Arguments and Outputs for details).
### E) Simulate changes in the spatial configuration of habitat in a landscape when patches are impermanent starting at QED
Calculate the **QED** of the number of patches present through time in that landscape, according to recipe E using the **make_reduced_Gmat( )** function. Then simply you may use the **simulate_impHabitat_BigN( )** function to run a simulation of patches coming and going through time, setting the parameters as desired (see #simulate_impHabitat_BigN for details).
#### E2) Calculate the average time until the habitat configuration will change from it’s current state
Use the **get_tau( )** function (see #get_tau in the functions section for details).
#### E3) Choose a new habitat configuration based on the probability the habitat will change from it’s current configuration to it
Supply the current habitat configuration and rates of change to the **transition( )** function. (See #transition for details).
### F) Calculate or estimate spatiotemporally weighted metrics of metapopulation size and persistence in a landscape of impermanent habitat (and other spatiotemporally weighted metrics)
First, if you haven’t got a **landscape** as out put by the **make_landscape( )** function do that now (see recipe A). Then consider the number of impermanent patches that in that landscape. If this number is small (<15-20), it is possible to calculate these metrics by supplying the output of recipe F1 to the **get_QEDMetapopMetrics( )** function with parameters set accordingly (see #get_QEDMetapopMetrics in functions for details), otherwise we can only estimate these metrics by supplying the output of recipe F3 which calculates weighted averages using only a subset of all the possible spatial configurations based on their probability at QED. **Calculate these once and save them since this is most often a computationally intensive task.**
#### F1) Calculate metrics of metapopulation size and persistence within all spatial configurations possible in a landscape of impermanent patches
 Construct every possible spatial configuration the habitat could be in by using the **make_allConfigs( )** function and supply the output to the **config_metrics_parallel( )** function with parameters set accordingly (see #config_metrics_parallel for details).
#### F2) Calculate metrics of metapopulation size and persistence within a subset of the spatial configurations possible in a landscape of impermanent patches
Construct a subset every possible spatial configuration the habitat could be in to supply to the **config_metrics_parallel( )** function with parameters set accordingly (see #config_metrics_parallel for details). 
#### F3) Estimate spatiotemporally weighted metrics of metapopulation size and persistence in a landscape of impermanent habitat (and other spatiotemporally weighted metrics)
To do this calculate the **QED** according to recipe E using the **make_reduced_Gmat( )** function. Supply this and how many spatial configurations of the habitat to use (**s.size**) in estimating these metrics (the closer this is to the number of all possible transient configurations, **2^n.patches-1**, the better but the longer and more memory this computation will take, easily exceeding the limits of even a decent server, so choose as large as is computationally feasible for how many of these calculations you want to perform and in what time frame) to the **subsample_reducedQED( )** function. Supply the output to recipe F2.
### 5) Simulate metapopulation dynamics in a landscape of impermanent patches
If you have not already done so make a **landscape** describing the spatial structure of habitat patches in the landscape (see recipe A). Then calculate the **QED** of the number of patches present through time in that landscape, according to recipe E using the **make_reduced_Gmat( )** function. Then simply you may use the **simulate_Metapop_impHabitat_BigN( )** function to run a simulation of a metapopulation in the landscape of impermanent habitat provided, setting the parameters as desired (see #simulate_Metapop_impHabitat_BigN for details).

___
## Glossary:
#### Landscape Parameters:
1) **n.patches:** takes an integer equal to the total number of patches to describe in a landscape
2) **landscape.type:** takes “line" or "plane" indicating whether patches are to be arranged either along a “line” or within a “plane” (1D or 2D landscape)
3) **landscape.size:** a numeric value indicating max y (and x coordinates if patches are arranged in 2D)
4) **p.dist:** takes “more clustered", "random", “more uniform”, or “uniform” indicating whether patches should be more clustered, randomly placed, more uniformly placed or uniformly placed with respect to each other within the landscape
5) **clust.its:** takes a positive integer equal to the number of times to iterate the algorithm increasing clustering or uniformity of the distribution of patches in the landscape if either the “more clustered” or “more uniform” option was chosen for the *p.dist* argument
6) **p.size.dist:** takes “lognormal", "random" or "uniform" indicating how patch sizes should be distributed within the landscape
7) **max.p.size:** numerical value indicating largest patch size
#### Impermanent Habitat Parameters:
2) **r.on:** a positive numeric value indicating the rate of patch recovery within a landscape
3) **r.off:** a positive numeric value indicating the rate of patch loss within a landscape
4) **s.size:** a positive integer indicating how many spatial configurations (excluding the absorbing no-habitat state) to take a subsample of from all of those possible (max 2^n.patches-1)
#### Metapopulation (Species Specific) Parameters:
7) **e.rate:** positive numeric value indicating the within patch population extinction rate
8) **c.rate:** positive numeric value indicating colonization rate when dispersers arrive at a patch
9) **delta:** *e.rate*/*c.rate*
10) **alpha:** a positive numeric value indicating 1/(avg. dispersal distance of a species)
11) **gamma:** a positive numeric value <1 indicating the strength to which upstream dispersal may be limited (0 no limitation), where 
12) **self.rec:** a positive numeric value indicating the extent to which disperses from a given patch may resettle on the same patch (0 makes resettlement impossible) iterations: number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
13) **iterations:** number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
#### Notable objects created and used by landscape functions:
8) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
9) **dist.mat:** a matrix of inter-patch distances ordered in rows and columns by their *patch.ID*. Note: upstream and downstream patches are not any more or less distant to each other based on stream order in this matrix *see disp.kernel* for how stream order can come into play.
#### Notable objects created and used by habitat CTMC functions:
10) **n.on:** a positive integer indicating the total number of patches present in a spatial configuration of habitat.
11) **config:** a binary vector indicating which patches (ordered by *patch.ID*) are present (*1*) or absent (*0*) within the landscape.
12) **configs** a binary sparse matrix where each column represents a unique patch which could be lost or recovered within a landscape and each row represents each unique spatial configuration that landscape could take on, where a 1 or 0 a patch as being either present or absent respectively.
13) **config_ID:** an positive integer identifying a unique spatial configuration of habitat *see Config_IDs document for details*
14) **G.mat:** the generator matrix of a CTMC describing random and independent patch recovery and loss from a landscape. Rows and columns correspond to either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches* x *2^n.patches*), or 2) how many patches are present in a spatial configuration (if *n.patches+1* x *n.patches+1*) *see Config_IDs document for details*.
15) **C.submat:** a submatrix of the generator matrix pertaining only to transitions amongst the transient (aka all excluding the absorbing no-habitat state) states of the habitat CTMC with rows and columns corresponding either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*.
16) **F.mat:** the Fundamental Matrix of the habitat CTMC with rows and columns corresponding  either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*.
17) **RMD:** the Ratio of Means Distribution of the habitat CTMC with rows and columns corresponding  either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*. The Ratio of Means Distribution gives the probabilities of states given the current state of a system and converges to the QED provided time to absorption (in this case all-off state) is long.
18) **QED:** a vector of positive numeric values indicating the probability of either 1) each unique spatial configurations of habitat ordered by *config_ID* (if length *2^n.patches-1*), or 2) how many patches are present in a spatial configuration or 3) each unique spatial configuration but indexed by how many patches are present in it (if length *n.patches*) at quasi-equilibrium *see Config_IDs document for details*. Which each function takes or outputs in the function descriptions. Note: the *get_QED( )* function returns not only this but also *rho1vrho2x2* in a list. By default only the first element of a list is used when passed to something that only takes one element.
19) **rho1vrho2x2:** a positive numeric value which when <1 indicates as per a rule of thumb ==cite== that the QED well describes the frequency of states as they are sampled regardless of the initial state.
20) **T.absorp:** a vector containing positive numeric values representing the average time to absorption (all-off state) from either 1) each unique spatial configurations of habitat ordered by *config_ID* (if length *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if length *n.patches*) *see Config_IDs document for details*.
21) **tau:** a positive numeric value of how much time will pass on average before the habitat configuration will change from it’s current state.
22) **configs.data:** a dataframe providing the configuration of the landscape (with which patches are on or off indicated by1's and 0's respectively) at each time point at which that configuration occurred
#### Notable objects created and used by the metapopulation (SRLM) functions
13) **disp.kernel:** a dataframe of values weighted relative to a species ability to move from one patch location to another (i.e. effective adjacency of patches) ordered in rows and columns by their *patch.ID* provided the *landscape* used to create this *disp.kernel* was also ordered by *patch.ID*.
14) **lambda.M:** a positive numeric value indicating the landscape capacity of a landscape (equals the metapopulation capacity lambda.M.delta if *delta=e.rate/c.rate=1*). The landscape is incapable of supporting a metapopulation with *delta=1* if *lambda.M<1*, but otherwise ensures persistence.
15) **lambda.M.delta:** a positive numeric value indicating the persistence capacity of a metapopulation in a landscape (equals *lambda.M* x *delta*). Also, represents the ==invasion== low density growth rate of the metapopulation for the SRLM. The landscape is incapable of supporting a metapopulation if *lambda.M<delta*, but otherwise ensures persistence.
16) **metapop.size:** a positive numeric value equal to the expected mean equilibrium occupancy of a landscape by a metapopulation (ranges from 0 indicating extinction to *n.patches* indicating the landscape is fully occupied).
17) **pstar:** a vector containing the equilibrium occupancy of each patch within a configuration (thus each subsequent column is ordered by patchID and contains the equilibrium occupancies of each patch).
#### Notable objects created and used by the metapopulation in impermanent habitat functions
14) **config.metrics:** a dataframe with columns for holding *lambda.M, lambda.M.delta, metapop.size, n.on, and pstar* calculated within either 1) a single spatial configuration of habitat (as done by the *get_MetapopMetrics()* function) providing 1 row of values for these metrics with that spatial configuration or 2) multiple spatial configurations of habitat (with rows ordered by the spatial configurations in which the metrics were calculated e.g. as supplied to the *get_config_metrics()* function).
15) **lm.QED.delta:** the geometric mean persistence capacity at QED.
16) **lm.QED:** the geometric mean landscape capacity at QED.
17) **geom.pstar.QED:** geometric mean equilibrium occupancy of the landscape by the metapopulation at QED.
18) **ar.pstar.QED:** arithmetic mean equilibrium occupancy of the landscape by the metapopulation at QED.
19) **exp.n.QED:** expected mean number of patches present at QED.
20) **pstar.configs:** a dataframe of the equilibrium occupancies of patches within each configuration supplied to this function (each column is ordered by patchID and contains the equilibrium occupancies of each patch within each configuration (ordered by row ordered by configurations supplied to the function).
21) **sim.data:** a dataframe providing the occupancy of each patch in the landscape over the increments of time on which the SRLM's system of ODE's is solved

___
## Functions:
*Note: The following information is also included in each function’s file respectively*
### Landscape Functions
### #make_landscape
#### DESCRIPTION
This function creates a data frame describing a number of patches in a landscape.
There are parameters to control:
1) how many patches there are to describe in a landscape
2) how the sizes (i.e. p.sizes or carrying capacities) of these patches vary within the landscape
3) how big to make the entire landscape 
4) whether patches are arranged along a 1D "line" or within a 2D "plane"
5) whether the patches tend to be "clustered" more together, versus "random"-ly or more "uniformal"-ly distributed throughout the landscape
#### USAGE
#### ARGUMENTS
1) **n.patches:** takes an integer equal to the total number of patches to describe in a landscape
2) **landscape.type:** takes “line" or "plane" indicating whether patches are to be arranged either along a “line” or within a “plane” (1D or 2D landscape)
3) **landscape.size:** a numeric value indicating max y (and x coordinates if patches are arranged in 2D)
4) **p.dist:** takes “more clustered", "random", “more uniform”, or “uniform” indicating whether patches should be more clustered, randomly placed, more uniformly placed or uniformly placed with respect to each other within the landscape
5) **clust.its:** takes a positive integer equal to the number of times to iterate the algorithm increasing clustering or uniformity of the distribution of patches in the landscape if either the “more clustered” or “more uniform” option was chosen for the *p.dist* argument
6) **p.size.dist:** takes “lognormal", "random" or "uniform" indicating how patch sizes should be distributed within the landscape
7) **max.p.size:** numerical value indicating largest patch size
#### DETAILS
#### VALUE
1) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
#### NOTE
### #get_distmat
#### DESCRIPTION
	This function extracts the inter-patch distance matrix from a *landscape* data frame. 
#### USAGE
#### ARGUMENTS
1) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
#### DETAILS
#### VALUE
1) **dist.mat:** a matrix of inter-patch distances ordered in rows and columns by their *patch.ID*. Note: upstream and downstream patches are not any more or less distant to each other based on stream order in this matrix *see disp.kernel* for how stream order can come into play.
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
#### NOTE
### Habitat CTMC Functions:
### #make_Gmat
#### DESCRIPTION
This function constructs the generator matrix (G.mat) for the CTMC, describing the transition rates from each habitat configuration (with a configuration ID corresponding to its row in G.mat) to another (with a configuration ID corresponding to its column in G.mat - see the Configuration IDs document for details). 
It does so by taking advantage of the block toeplitz structure of the generator matrix to save on computation and memory - see the Configuration IDs document for details
#### USAGE
#### ARGUMENTS
1) **n.patches:** takes an integer equal to the total number of patches to describe in a landscape
2) **r.on:** a positive numeric value indicating the rate of patch recovery within a landscape
3) **r.off:** a positive numeric value indicating the rate of patch loss within a landscape
#### DETAILS
#### VALUE
1) **G.mat:** the generator matrix of a CTMC describing random and independent patch recovery and loss from a landscape. Rows and columns correspond to each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches* x *2^n.patches*) *see Config_IDs document for details*.
#### REQUIRED PACKAGES
* **Matrix:** to get sparse matrix functions
* **gdata:** to get upper and lower triangles of matrices
* **lava:** to get the anti-diagonal of a matrix using revdiag
#### REQUIRED FUNCTIONS
none
#### NOTE
### #make_reduced_Gmat
#### DESCRIPTION
This function constructs the generator matrix (*G.mat*) for the reduced CTMC, describing the transition rates from habitat configurations with some number of patches on (corresponding to the row in G.mat) to some those with some other number of patches on (corresponding to the column in G.mat - see the Configuration IDs document for details). 
It does so by taking advantage of the toeplitz structure of the generator matrix to save on computation and memory - see the Configuration IDs document for details
**NOTE: the following warning messages are suppressed in the function code because they inconsequential and to be ignored: 
Warning messages:
1: In x[upper.tri(x, diag = diag)] <- value :number of items to replace is not a multiple of replacement length
2: In ret[rev(lower.tri(x, diag = diag))] <- value :number of items to replace is not a multiple of replacement length**
#### USAGE
#### ARGUMENTS
1) **n.patches:** takes an integer equal to the total number of patches to describe in a landscape
2) **r.on:** a positive numeric value indicating the rate of patch recovery within a landscape
3) **r.off:** a positive numeric value indicating the rate of patch loss within a landscape
#### DETAILS
#### VALUE
1) **G.mat:** the generator matrix of a CTMC describing random and independent patch recovery and loss from a landscape. Rows and columns correspond to how many patches are present in a spatial configuration (if *n.patches+1* x *n.patches+1*) *see Config_IDs document for details*.
#### REQUIRED PACKAGES
* **Matrix:** to get sparse matrix functions
* **gdata:** to get upper and lower triangles of matrices
* **lava:** to get the anti-diagonal of a matrix using revdiag
#### REQUIRED FUNCTIONS
none
#### NOTE
### #Gmat_to_Csubmat
#### DESCRIPTION
This function simply removes transitions into and out of the absorbing (all-off) state the CTMC's generator matrix (G.mat) to provide the submatrix (C.mat), describing only the transition rates among the transient (non-absorbing states). Works the same regardless of whether G.mat is describing transitions among states corresponding to the number of patches on or among specific configurations of patches on and off since the first row and column always corresponds to the absorbing (all-off) state.
#### USAGE
#### ARGUMENTS
1) **G.mat:** the generator matrix of a CTMC describing random and independent patch recovery and loss from a landscape. Rows and columns correspond to either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches* x *2^n.patches*), or 2) how many patches are present in a spatial configuration (if *n.patches+1* x *n.patches+1*) *see Config_IDs document for details*.
#### DETAILS
#### VALUE
1) **C.submat:** a submatrix of the generator matrix pertaining only to transitions amongst the transient (aka all excluding the absorbing no-habitat state) states of the habitat CTMC with rows and columns corresponding either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*.
#### REQUIRED PACKAGES:
none
#### REQUIRED FUNCTIONS:
none
#### NOTE
### #get_QED
#### DESCRIPTION
This function calculates the quasi-equilibrium distribution (Q.E.D. or QED) which provides the frequency of each of the transient states while not in the absorbing (all-off) state, regardless of the initial state provided convergence to the Q.E.D. is fast relative to the absorbing (all-off state). It also calculates a rule of thumb to check if this is the case.
**Note: By default only the first element of a list is used when passed to something that only takes one element.**
#### USAGE
#### ARGUMENTS
1) **C.submat:** a submatrix of the generator matrix pertaining only to transitions amongst the transient (aka all excluding the absorbing no-habitat state) states of the habitat CTMC with rows and columns corresponding either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*.
#### DETAILS
#### VALUE
1) a list with:
   1) **QED:** a vector of positive numeric values indicating the probability of either 1) each unique spatial configurations of habitat ordered by *config_ID* (if length *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if length *n.patches*) at quasi-equilibrium *see Config_IDs document for details*. 
   2) **rho1vrho2x2:** a positive numeric value which when <1 indicates as per a rule of thumb ==cite== that the QED well describes the frequency of states as they are sampled regardless of the initial state.
#### REQUIRED PACKAGES
* **Matrix:** for working with sparse matrices
#### REQUIRED FUNCTIONS
#### NOTE
### #Csubmat_to_Fmat
#### DESCRIPTION
This function uses the generator matrix (*G.mat*) for the CTMC to calculate the Fundamental Matrix (*F.mat*), mathematically defined as the negative inverse of *C.mat*
#### USAGE
#### ARGUMENTS
1) **C.submat:** a submatrix of the generator matrix pertaining only to transitions amongst the transient (aka all excluding the absorbing no-habitat state) states of the habitat CTMC with rows and columns corresponding either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details* 
#### DETAILS
#### VALUE
3) **F.mat:** the Fundamental Matrix of the habitat CTMC with rows and columns corresponding  either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*
#### REQUIRED PACKAGES
* **Matrix:** for working with sparse matrices
#### REQUIRED FUNCTIONS
none
#### NOTE
### #Fmat_to_RMD
#### DESCRIPTION
This function uses the fundamental matrix (*F.mat*) for the CTMC to calculate the ratio of means distribution (RMD), aka the frequency of each state (given by each column) provided the system began in a given state (given by each row) before absorption (ending in  the all-off state). The ratio of means distribution converges to the Q.E.D. provided time to absorption is long enough and thus, can be used as a way to check the CTMC converges to the Q.E.D. sufficiently quickly relative to the absorbing state, such that the Q.E.D. describes the frequency of each state regardless of what state the CTMC may have began in.
#### USAGE
#### ARGUMENTS
1) **F.mat:** the Fundamental Matrix of the habitat CTMC with rows and columns corresponding  either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*
#### DETAILS
#### VALUE
1) **RMD:** the Ratio of Means Distribution of the habitat CTMC with rows and columns corresponding  either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*. The Ratio of Means Distribution gives the probabilities of states given the current state of a system and converges to the QED provided time to absorption (in this case all-off state) is long.
#### REQUIRED PACKAGES
none
#### NOTE
### #get_Tabsorp
#### DESCRIPTION
This function calculates the average time to absorption (all-off state) from each state (organized by row) using *C.submat*. Mathematically, this is given by summing across rows of the fundamental matrix and multiplying by the trace of (I-*C.submat*) (i.e. the rate of leaving states, where I is an identity matrix) **Note: This requires calculating the fundamental matrix which for large numbers of patches is a computationally costly endeavor (should you wish to calculate this and the ratio of means distribution many times over it would be better to first calculate *F.mat* outside this function and then use it as an argument to run both the *Fmat_to_RMD* function and this function (adding it as an argument and commenting out step 1). I have not bothered to do this here as I only intend to calculate the ratio of means distribution for a few test cases just for checking coding and convergence to Q.E.D.**
#### USAGE
#### ARGUMENTS
1) **C.submat:** a submatrix of the generator matrix pertaining only to transitions amongst the transient (aka all excluding the absorbing no-habitat state) states of the habitat CTMC with rows and columns corresponding either 1) each unique spatial configurations of habitat ordered by *config_ID* (if *2^n.patches-1* x *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if *n.patches* x *n.patches*) *see Config_IDs document for details*
#### DETAILS
#### VALUE
1) **T.absorp:** a vector containing positive numeric values representing the average time to absorption (all-off state) from either 1) each unique spatial configurations of habitat ordered by *config_ID* (if length *2^n.patches-1*), or 2) how many patches are present in a spatial configuration (if length *n.patches*) *see Config_IDs document for details*.
#### REQUIRED PACKAGES
* **Matrix:** for working with sparse matrices
#### REQUIRED FUNCTIONS
none
#### NOTE
### #get_tau
#### DESCRIPTION
This function draws a value for when the system will leave a given habitat configuration (i.e. how long until the system will on average jump to a new configuration from the current configuration)
#### USAGE
#### ARGUMENTS
1) **config:** a binary vector indicating which patches (ordered by *patch.ID*) are present (*1*) or absent (*0*) within the landscape.
2) **r.on:** a positive numeric value indicating the rate of patch recovery within a landscape
3) **r.off:** a positive numeric value indicating the rate of patch loss within a landscape
#### DETAILS
#### VALUE
1) **tau:** a positive numeric value of how much time will pass on average before the habitat configuration will change from it’s current state.
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
#### NOTE
### #transition
#### DESCRIPTION
This function makes a transition from one habitat configuration to another based on the probability of transitioning to that configuration from the current configuration given by the generator matrix
#### USAGE
#### ARGUMENTS
1) **config:** a binary vector indicating which patches (ordered by *patch.ID*) are present (*1*) or absent (*0*) within the landscape.
2) **r.on:** a positive numeric value indicating the rate of patch recovery within a landscape
3) **r.off:** a positive numeric value indicating the rate of patch loss within a landscape
#### DETAILS
#### VALUE
1) **config:** a binary vector indicating which patches (ordered by *patch.ID*) are present (*1*) or absent (*0*) within the landscape.
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
#### NOTE
### #IDs_to_Configs
#### DESCRIPTION
This function convert identifying numbers for configurations into the configurations (binary vectors indicating which patches are on and off) they represent. It uses an algorithm detailed in the Configuration IDs document.
#### USAGE
#### ARGUMENTS
1) **config_ID:** an positive integer identifying a unique spatial configuration of habitat *see Config_IDs document for details*
2) **n.patches:** takes an integer equal to the total number of patches to describe in a landscape
#### DETAILS
#### VALUE
1) **config:** a binary vector indicating which patches (ordered by *patch.ID*) are present (*1*) or absent (*0*) within the landscape.
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
#### NOTE
### #make_allConfigs
#### DESCRIPTION
This function constructs a sparse matrix where each row represents each unique configuration a landscape with a given number of patches could take on, where a 1 or 0 in each column indicates whether a patch identified by that column being either present or absent within a given landscape of those patches. 
The row order of configurations corresponds to their row and column order in the generator matrix (*G_mat*) and likewise, the construction of habitat configurations takes advantage of block structure by which the generator matrix grows to save on computation and memory - **see Config_IDs document for details.**
**While this function is extremely computationally and memory efficient, the task completed by this function becomes impossible for landscapes with large numbers of patches (n.patches) because the number of configurations possible grows to 2^n.patches and thus becomes huge very quickly!**
#### USAGE
#### ARGUMENTS:
2) **n.patches:** takes an integer equal to the total number of patches to describe in a landscape
#### DETAILS
#### VALUE
1) **configs:** a binary sparse matrix where each column represents a unique patch which could be lost or recovered within a landscape and each row represents each unique spatial configuration that landscape could take on, where a 1 or 0 a patch as being either present or absent respectively.
#### REQUIRED PACKAGES
* **Matrix:** to get sparse matrices and vectors
#### REQUIRED FUNCTIONS
none
#### NOTE
### #reduced_to_full_QED
#### DESCRIPTION
Because patches randomly and independently are gained and lost, the probabilities of habitat configurations with the same number of habitat patches present are equal at QED. Therefore, the QED of the reduced CTMC describing only the number of habitat patches on converts to the QED of the full CTMC describing every possible habitat configuration by dividing those probabilities by the number of those habitat configurations with the same number of habitat patches on.
This function uses this property to covert the QED probabilities of states with a given number of patches on to the QED probabilities of each habitat configuration with that number of habitat patches on. 
#### USAGE
#### ARGUMENTS
10) **QED:** a vector of positive numeric values indicating the probability of how many patches are present in a spatial configuration (length *n.patches*) at quasi-equilibrium *see Config_IDs document for details*.
11) **n.patches:** takes an integer equal to the total number of patches to describe in a landscape
#### DETAILS
#### VALUE
10) **QED:** a vector of positive numeric values indicating the probability of each unique spatial configuration but indexed by how many patches are present in it (length *n.patches*) at quasi-equilibrium *see Config_IDs document for details*.
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
#### NOTE
### #subsample_reducedQED
#### DESCRIPTION
This function constructs the configurations (binary vectors indicating which patches are on and off) for just a subset of all the configurations possible for a given number of patches chosen according to the probability with which these configurations occur at Q.E.D. It uses an algorithm detailed in the Configuration IDs document.
#### USAGE
#### ARGUMENTS
10) **QED:** a vector of positive numeric values indicating the probability of how many patches are present in a spatial configuration (length *n.patches*) at quasi-equilibrium *see Config_IDs document for details*.
11) **s.size:** a positive integer indicating how many spatial configurations (excluding the absorbing no-habitat state) to take a subsample of from all of those possible (max 2^n.patches-1)
#### DETAILS
#### VALUE
1) **configs** a binary sparse matrix where each column represents a unique patch which could be lost or recovered within a landscape and each row represents each unique spatial configuration that landscape could take on, where a 1 or 0 a patch as being either present or absent respectively.
#### REQUIRED PACKAGES
* **tidyverse:** to get filter function
* **oddeven:** to get even() function
#### REQUIRED FUNCTIONS
* IDsToConfigs.r
* reduced_to_full_QED
#### NOTE
### #simulate_impHabitat_BigN
#### DESCRIPTION
This function simulates the dynamics of a landscape of impermanent habitat patches gained and lost according to a CTMC describing how the landscape's configurations of available habitat patches change over time. The code to run these simulations has been extensively optimized for speed and memory usage and can work even when there is a large number of patches in a landscape (e.g. *n.patches > 20*), but nonetheless computational limits may be reached if the number of patches is too large. 
#### USAGE
#### ARGUMENTS
1) **limit:** takes a positive integer indicating either the the maximum number of transitions or duration of time the simulation should be run for
2) **limit.type:** takes “trans” or “time” for whether *limit* indicates the maximum number of transitions or duration of time the simulation should be run for respectively, the default is “time”
3) **r.on:** a positive numeric value indicating the rate of patch recovery within a landscape
4) **r.off:** a positive numeric value indicating the rate of patch loss within a landscape
5) **QED:** a vector of positive numeric values indicating the probability of how many patches are present in a spatial configuration (length *n.patches*) at quasi-equilibrium *see Config_IDs document for details*.
#### DETAILS
#### VALUE 
1) **configs.data:** a dataframe providing the configuration of the landscape (with which patches are on or off indicated by1's and 0's respectively) at each time point at which that configuration occurred
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
* transition
* subsample_reducedQED
* reduced_to_full_QED
* IDs_to_configs: IDsToConfigs
* get_tau
#### NOTE
### Metapopulation SRLM Functions:
### #get_dispkernel
#### DESCRIPTION
This function calculates a dispersal kernel as it might be specified for a species' metapopulation
#### USAGE
#### ARGUMENTS
1) **dist.mat:** a matrix of inter-patch distances ordered in rows and columns by their *patch.ID*. Note: upstream and downstream patches are not any more or less distant to each other based on stream order in this matrix *see disp.kernel* for how stream order can come into play.
2) **alpha:** a positive numeric value indicating 1/(avg. dispersal distance of a species)
3) **gamma:** a positive numeric value <1 indicating the strength to which upstream dispersal may be limited (0 no limitation), where 
4) **self.rec:** a positive numeric value indicating the extent to which disperses from a given patch may resettle on the same patch (0 makes resettlement impossible) iterations: number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
#### DETAILS
#### VALUE
13) **disp.kernel:** a dataframe of values weighted relative to a species ability to move from one patch location to another (i.e. effective adjacency of patches) ordered in rows and columns by their *patch.ID* provided the *landscape* used to create this *disp.kernel* was also ordered by *patch.ID*.
#### REQUIRED PACKAGES:
none
#### REQUIRED FUNCTIONS:
none
#### NOTE
### #get_MetapopMetrics
#### ==DESCRIPTION==
#### USAGE
#### ARGUMENTS
1) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
2) **e.rate:** positive numeric value indicating the within patch population extinction rate
3) **c.rate:** positive numeric value indicating colonization rate when dispersers arrive at a patch
4) **delta:** *e.rate*/*c.rate*
5) **alpha:** a positive numeric value indicating 1/(avg. dispersal distance of a species)
6) **gamma:** a positive numeric value <1 indicating the strength to which upstream dispersal may be limited (0 no limitation), where 
7) **self.rec:** a positive numeric value indicating the extent to which disperses from a given patch may resettle on the same patch (0 makes resettlement impossible) iterations: number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
8) **iterations:** number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
#### DETAILS
#### VALUE
1) a list of objects containing:
   1) **lambda.M:** a positive numeric value indicating the landscape capacity of a landscape (equals the metapopulation capacity lambda.M.delta if *delta=e.rate/c.rate=1*). The landscape is incapable of supporting a metapopulation with *delta=1* if *lambda.M<1*, but otherwise ensures persistence.
   2) **lambda.M.delta:** a positive numeric value indicating the persistence capacity of a metapopulation in a landscape (equals *lambda.M* x *delta*). Also, represents the ==invasion== low density growth rate of the metapopulation for the SRLM. The landscape is incapable of supporting a metapopulation if *lambda.M<delta*, but otherwise ensures persistence.
   3) **metapop.size:** a positive numeric value equal to the expected mean equilibrium occupancy of a landscape by a metapopulation (ranges from 0 indicating extinction to *n.patches* indicating the landscape is fully occupied).
   4) **pstar:** a vector containing the equilibrium occupancy of each patch within a configuration (thus each subsequent column is ordered by patchID and contains the equilibrium occupancies of each patch).
#### REQUIRED PACKAGES
none 
#### REQUIRED FUNCTIONS
* get_lambdaM
* get_pstar
* get_distmat
* get_dispkernel
#### NOTE
### #get_lambdaM
#### DESCRIPTION
This function calculates the landscape or metapopulation capacity (lambda.M, whether e.rate/c.rate=1 or not respectively) within a given configuration of habitat
#### USAGE
#### ARGUMENTS
1) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
2) **e.rate:** positive numeric value indicating the within patch population extinction rate
3) **c.rate:** positive numeric value indicating colonization rate when dispersers arrive at a patch
4) **delta:** *e.rate*/*c.rate*
5) **disp.kernel:** a dataframe of values weighted relative to a species ability to move from one patch location to another (i.e. effective adjacency of patches) ordered in rows and columns by their *patch.ID* provided the *landscape* used to create this *disp.kernel* was also ordered by *patch.ID*.
#### DETAILS
#### VALUE
1) ****lambda.M:** a positive numeric value indicating the landscape capacity of a landscape (equals the metapopulation capacity lambda.M.delta if *delta=e.rate/c.rate=1*). The landscape is incapable of supporting a metapopulation with *delta=1* if *lambda.M<1*, but otherwise ensures persistence.
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
#### NOTE
### #get_pstar
#### DESCRIPTION
This function calculates the equilibrium occupancy (probabilities with which patches occupied at equilibrium by a metapopulation) within a given configuration of habitat
#### USAGE
#### ARGUMENTS
1) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
2) **e.rate:** positive numeric value indicating the within patch population extinction rate
3) **c.rate:** positive numeric value indicating colonization rate when dispersers arrive at a patch
4) **delta:** *e.rate*/*c.rate*
5) **disp.kernel:** a dataframe of values weighted relative to a species ability to move from one patch location to another (i.e. effective adjacency of patches) ordered in rows and columns by their *patch.ID* provided the *landscape* used to create this *disp.kernel* was also ordered by *patch.ID*.
6) **iterations:** number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
#### DETAILS
#### VALUE
1) **pstar:** a vector containing the equilibrium occupancy of each patch within a configuration (thus each subsequent column is ordered by patchID and contains the equilibrium occupancies of each patch).
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
#### NOTE
### #name_SRLMODEparams
#### DESCRIPTION
This function both creates 1) a named list of the parameters to be supplied to the SRLMODE function because R’s ode solvers require all parameters to be supplied in a single named list of with unique names for every parameter use in an equation (aka if a parameter is a vector or a matrix every single element of that vector or matrix must have a unique parameter name! how annoying!) and 2) the names it should to use for the interpatch distances if these do not already exist and haven't been supplied to this function.
#### USAGE
#### ARGUMENTS
1) **parameter.values:** a concatenated list of containing all the objects defining all the parameter values to be used to construct the system of differential equations for the SRLM: n.patches, c.rate, self.rec, extinction.rates, x.coord, y.coord, areas, config, f.vals
#### DETAILS
#### VALUE
1) a list containing:
   1) **parameters:** a named list where each entry contains the unique name of a parameter and its value, for every single value needed to construct the system of differential equations for the SRLM 
   2) **f.names:** a vector of character values providing a "name" for each interpatch distance
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS 
name_Mij
#### NOTE
### #SRLMODE
#### DESCRIPTION
This function simply defines the system of differential equations describing the SRLM as they must be set up for use by one of R's ode solvers.
#### USAGE
#### ARGUMENTS
1) **t:** a time sequence for the equations to be solved over
2) **p:** occupancy of patches as described by the equations set up by this function
3) **parameters:** a vector of parameters used to construct the equations (must be a 
named list with unique names for every parameter, including all elements of 
any vectors or matrices)
#### USAGE
### #at_equilibrium_rootfunc
Provides a stopping condition to solving the SRLM when the change in patch occupancy is less than 1e-4 and therefore the SRLM is essentially at equilibrium
#### USAGE
#### ARGUMENTS
1) **t:** a time sequence for the equations to be solved over
2) **p:** occupancy of patches as described by the equations set up by this function
3) **parameters:** a vector of parameters used to construct the equations (must be a 
named list with unique names for every parameter, including all elements of 
any vectors or matrices)
#### NOTE
### Metapopulations in Impermanent Habitat Functions:
### #get_config_metrics
#### DESCRIPTION
This function calculates the landscape capacity (*lambda.M*), persistence capacity (*lambda.M.delta*), equilibrium occupancy of each patch (*p.star*), total equilibrium occupancy across patches in the landscape (*metapop.size*), and the number of patches on (n.on) within a single habitat configuration in a dynamic landscape
#### USAGE
#### ARGUMENTS
1) **config:** a binary vector indicating which patches (ordered by *patch.ID*) are present (*1*) or absent (*0*) within the landscape.
2) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
3) **e.rate:** positive numeric value indicating the within patch population extinction rate
4) **c.rate:** positive numeric value indicating colonization rate when dispersers arrive at a patch
5) **delta:** *e.rate*/*c.rate*
6) **alpha:** a positive numeric value indicating 1/(avg. dispersal distance of a species)
7) **gamma:** a positive numeric value <1 indicating the strength to which upstream dispersal may be limited (0 no limitation), where 
8) **self.rec:** a positive numeric value indicating the extent to which disperses from a given patch may resettle on the same patch (0 makes resettlement impossible) iterations: number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
#### DETAILS
#### VALUE
10) **config.metrics:** a dataframe with columns for holding *lambda.M, lambda.M.delta, metapop.size, n.on, and pstar* calculated within  a single spatial configuration of habitat.
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
* get_distmat
* get_dispkernel
* get_lambdaM
* get_pstar
#### NOTE
### #config_metrics_parallel
#### DESCRIPTION
This function calculates the landscape capacity (*lambda.M*), persistence capacity (*lambda.M.delta*), equilibrium occupancy of each patch (*p.star*), total equilibrium occupancy across patches in the landscape (*metapop.size*), and the number of patches on (n.on) and places these in a dataframe for some subset of all possible habitat configurations of a dynamic landscape assigned to be performed by a single core, when calculating these metrics over some number of cores (*N*) in parallel.
**Note: To use this function you must update the file path for where the cores should source the required functions on the server you are using.**
#### USAGE
#### ARGUMENTS
1) **I:** an index indicating which core on which to perform the calculations
2) **N:** the number of cores over which the calculations are being performed in parallel over
3) **configs:** a binary sparse matrix where each column represents a unique patch which could be lost or recovered within a landscape and each row represents each unique spatial configuration that landscape could take on, where a 1 or 0 a patch as being either present or absent respectively.
4) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
5) **iterations:** number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
#### DETAILS
#### VALUE 
14) **config.metrics:** a dataframe with columns for holding *lambda.M, lambda.M.delta, metapop.size, n.on, and pstar* calculated within multiple spatial configurations of habitat (with rows ordered by the spatial configurations in which the metrics were calculated e.g. as supplied to the *get_config_metrics()* function).
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
* get_distmat
* get_dispkernel
* get_lambdaM
* get_pstar
#### NOTE
### #get_QEDMetapopMetrics
#### DESCRIPTION
This function computes expectations of metapopulation size and persistence, and number of habitat patches at QED by taking the landscape capacity (lambda.M), persistence capcity (lambda.M.delta), equilibrium occupancy (metapop.size), and number of habitat patches present within habitat configurations and computing weighted averages based on the probabilities with which these habitat configurations occur at QED.
#### USAGE
#### ARGUMENTS
14) **config.metrics:** a dataframe with columns for holding *lambda.M, lambda.M.delta, metapop.size, n.on, and pstar* calculated within spatial configurations of habitat (with rows ordered by the spatial configurations in which the metrics were calculated e.g. as supplied to the *get_config_metrics()* function). **If these metrics weren’t calculated for every possible spatial configuration then this function can only provide estimates of it’s output (more spatial configurations used = better estimates)**
#### DETAILS
#### VALUE
1) a list of objects containing:
   1) **lm.QED.delta:** the geometric mean persistence capacity at QED.
   2) **lm.QED:** the geometric mean landscape capacity at QED.
   3) **geom.pstar.QED:** geometric mean equilibrium occupancy of the landscape by the metapopulation at QED.
   4) **ar.pstar.QED:** arithmetic mean equilibrium occupancy of the landscape by the metapopulation at QED.
   5) **exp.n.QED:** expected mean number of patches present at QED.
   6) **pstar.configs:** a dataframe of the equilibrium occupancies of patches within each configuration supplied to this function (each column is ordered by patchID and contains the equilibrium occupancies of each patch within each configuration (ordered by row ordered by configurations supplied to the function).
#### REQUIRED PACKAGES
none 
#### REQUIRED FUNCTIONS
none
#### NOTE
### #simulate_Metapop_impHabitat_BigN
#### DESCRIPTION
This function simulates the dynamics of a metapopulation in a landscape of impermanent habitat patches gained and lost according to a CTMC describing how the landscape's configurations of available habitat patches change over time. The code to run these simulations has been extensively optimized for speed and memory usage and can work even when there is a large number of patches in a landscape (e.g. *n.patches > 20*), but nonetheless computational limits may be reached if the number of patches is too large. 
#### USAGE
#### ARGUMENTS
1) **limit:** takes a positive integer indicating either the the maximum number of transitions or duration of time the simulation should be run for
2) **limit.type:** takes “trans” or “time” for whether *limit* indicates the maximum number of transitions or duration of time the simulation should be run for respectively, the default is “time”
3) **landscape:** a data frame containing patch.ID, p.sizes, coordinates, and a distance.matrix with columns formatted as follows
   - **patch.ID:** is a unique integer ID number for each patch. If there is a stream order to patches, these run from downstream to upstream and it is important to not change the order of rows in the landscape dataframe when calculating a dispersal kernel for a species in the landscape.
   - **p.sizes:** provides the size of each patch
   - **y.coord:** provides vertical (e.g. latitudinal) patch locations on a "map"
   - **x.coord:** provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
   - subsequent columns contain inter-patch distances forming a distance matrix (see *dist.mat* for details)
4) **e.rate:** positive numeric value indicating the within patch population extinction rate
5) **c.rate:** positive numeric value indicating colonization rate when dispersers arrive at a patch
6) **delta:** *e.rate*/*c.rate*
7) **alpha:** a positive numeric value indicating 1/(avg. dispersal distance of a species)
8) **gamma:** a positive numeric value <1 indicating the strength to which upstream dispersal may be limited (0 no limitation), where 
9) **self.rec:** a positive numeric value indicating the extent to which disperses from a given patch may resettle on the same patch (0 makes resettlement impossible) iterations: number of iterations for which the iterative function f for finding pstar is iterated -greater iterations = greater accuracy, less iterations = lower accuracy
10) **r.on:** a positive numeric value indicating the rate of patch recovery within a landscape
11) **r.off:** a positive numeric value indicating the rate of patch loss within a landscape
12) **QED:** a vector of positive numeric values indicating the probability of how many patches are present in a spatial configuration (length *n.patches*) at quasi-equilibrium *see Config_IDs document for details*.
#### DETAILS
#### VALUE 
1) a list of
   1) **sim.data:** a dataframe providing the occupancy of each patch in the landscape over the increments of time on which the SRLM's system of ODE's is solved
   2) **configs.data:** a dataframe providing the configuration of the landscape (with which patches are on or off indicated by1's and 0's respectively) at each time point at which that configuration occurred
#### REQUIRED PACKAGES
* deSolve: for ode solver
#### REQUIRED FUNCTIONS
* transition
* SRLMODE
* at_equilibrium_rootfunc
* get_distmat
* get_dispkernel
* get_pstar
* name_SRLMODEparams
* subsample_reducedQED
* reduced_to_full_QED
* IDs_to_configs: IDsToConfigs
* get_tau
* get_Mij
#### NOTE
==ADD PARALLEL SIMS FUNCTION==
### Additional Functions made simply because they are useful:
### #split_jobs
#### DESCRIPTION
Assuming you have an indexed set of jobs (numbered 1-J) to be run over only N computational cores, this function splits that set into an indexed subset of those jobs to be run on a core. e.g. if I have 10 jobs numbered 1-10 to run in parallel and only 3 computational cores available this function will give jobs 1, 2, 3, and 4 to core 1, jobs 5, 6, and 7 to core 2, and jobs 8, 9, 10 to core 3. 
N is the maximum number of core you would like the jobs split over. If the number of jobs (J) is less than N, unnecessary cores will remain unused. If I had 2 jobs, core gets 1 job, core 2 gets job 2, and core 3 is unused.
#### USAGE
#### ARGUMENTS
1) **N:** the maximum number of cores to divide the jobs over
2) **J:** the number of jobs to divide over them
3) **i:** the core index (i.e. index indicating which core to use)
#### DETAILS
#### VALUE
1) **jobs.for.i:** a vector providing an index of the jobs to be run on core i
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS
none
### #name_Mij
#### DESCRIPTION
This function creates unique names for each element of a matrix.
#### USAGE
#### ARGUMENTS
1) **mat:** a matrix
2) **mat.name:** a character string providing a name for the matrix (each element will be named by this followed by it's row and column)
#### DETAILS
#### VALUE
1) **m.names:** a matrix of character strings providing names for each element of the matrix where each element will be named by this followed by it's row and column.
#### REQUIRED PACKAGES
none
#### REQUIRED FUNCTIONS 
none
#### NOTE

___

## Scripts:
### Script for calculating or estimating QEDMetapopMetrics with a subset of configs (large N):
Source this script to calculate the persistence capacity (lambda.M.delta), landscape capacity (lambda.M), geometric mean equilibrium occupancy of the landscape by the metapopulation at QED (geom.pstar.QED), arithmetic mean equilibrium occupancy of the landscape by the metapopulation at QED (ar.pstar.QED), expected mean number of patches present at QED (exp.n.QED), and a dataframe of the equilibrium occupancies of patches within each configuration used in the calculation of these metrics. It does so for a metapopulation (as specified by the metapop parameters *see below) in a dynamic landscape of impermanent patches with a given QED of habitat configurations but only for a given set of those habitat configurations (with how many to use specified by the parameter s.size). This set may include up to all possible 2^n.patches configurations (when computationally feasible) or a subset of these configurations chosen by their probability of occurrence at QED to provide estimates of these metrics in landscapes containing larger numbers of patches (when using all 2^n.patches configurations is computationally infeasible). To increase computational efficiency within habitat configuration calculations are performed in parallel over N cores (where N is the number of cores available). However, computation of these metrics may still be limited by the amount of memory available even using a relatively small number of habitat configurations possible in landscapes containing larger numbers of patches. This is because the algorithm used to convert ID #'s associated with each configuration possible in a landscape still requires storage of vectors at least as long as the number of configurations possible with the same number of patches. This is still much less than 2^n.patches but these still become big quickly as the number of patches in landscapes increases. e.g. Using a server with ~ 50 cores and lots of memory it is very feasible to perform these computations with landscapes of ~30 patches, but 50 patches becomes computationally impossible. Using a larger percentage of all the configurations possible in a landscape increases the accuracy of these metrics but increases the time required to compute these metrics based on the number of cores available to you. e.g. If you have ~50 cores available for these calculations then calculations for ~50x more configurations can be used and computed in the same time it would take to compute these metrics using the same number of configurations on 1 core. #ConfigIDs_doc

### Script for calculating QEDMetapopMetrics using all possible configs (small N):
**WARNING: SLOW & MEMORY INTENSIVE! USE ONLY FOR SMALL LANDSCAPES (n.patches~<=20)**
Source this script to calculate the persistence capacity (lambda.M.delta), landscape capacity (lambda.M), geometric mean equilibrium occupancy of the landscape by the metapopulation at QED (geom.pstar.QED), arithmetic mean equilibrium occupancy of the landscape by the metapopulation at QED (ar.pstar.QED), expected mean number of patches present at QED (exp.n.QED), and a dataframe of the equilibrium occupancies of patches within every habitat configuration possible.
Because the number of habitat configurations possible becomes huge (size 2^n.patches) for large numbers of patches (n.patches), this script can only be run for landscapes with a small number of patches.

==ADD PARALLEL SIMS SCRIPT==



