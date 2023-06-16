# msam_australian_birds

This repository hosts data and R code to reproduce the Multi-Species Abundance Model (MSAM) from García-Navas V., López-Poveda G., Bliard L., Christidis L., Ozgul A. (2023). No country for small birds: positive association among medium-sized, aggressive species in Australian bird communities.


## GENERAL INFORMATION

1. Title: "Data and code to reproduce the Multi-Species Abundance Model (MSAM) from: No country for small birds, positive association among medium-sized, aggressive species in Australian bird communities"

2. Author information:
       
       Name: Louis Bliard
		   Institution: University of Zurich
		   Email: Louis.bliard@evobio.eu

3. Date of data collection (single date, range, approximate date): 1998-2020

4. Geographic location of data collection: south-east Queensland

5. Information about funding sources that supported the collection of the data: NA


## SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data:

2. Links to publications that cite or use the data: García-Navas V, López-Povedas G, Bliard L, Chrsitidis L, Ozgul A. 2023. No country for small birds: positive association among medium-sized, aggressive species in Australian bird communities.

3. Links to other publicly accessible locations of the data: Data associated to "García-Navas V, López-Povedas G, Bliard L, Chrsitidis L, Ozgul A. 2023. No country for small birds: positive association among medium-sized, aggressive species in Australian bird communities". Figshare. https://doi.org/10.6084/m9.figshare.20463717.v1

4. Links/relationships to ancillary data sets: NA

5. Was data derived from another source? yes
	A. If yes, list source(s): Australian Bird Atlas https://birdata.birdlife.org.au


## METHODOLOGICAL INFORMATION


1. Methods for processing the data: All data were processed in R.

2. Software-specific information needed to interpret the data:
R statistical software, version 4.1.2. 
JAGS https://mcmc-jags.sourceforge.io/


## DATA & FILE OVERVIEW

1. File List: 

- `data_and_output.RData` R data object containing the bird surveys data needed to run the MSAM model, and the output of the model.

- `msam_analysis.R` Code containing the JAGS model to run the MSAM analysis, which is used to compute estimated true abundance for species at different sites.

2. Relationship between files, if important: 
The RData file "data_and_output.RData" is needed to run the MSAM analysis with the "msam_analysis.R" file. 


### DATA-SPECIFIC INFORMATION FOR: `data_and_output.RData`

1. Number of variables: NA

2. Number of cases/rows: NA

3. Object List: 
- f1: vector of length 52 containing the season of first survey for each site.
- f2: vector of length 52 containing the season of last survey for each site.
- Nmat: 4 dimensional array of dimensions 52 (number of sites) * 17 (number of seasons) * 134 (number of species) * 600 (number of saved iterations). This is the output of the MSAM model, and represents the estimated abundance of each species at each site for each season.
- season_vec: vector of length 17 corresponding to each season (starting from season 1998/1999).
- site_vec: vector of length 52 corresponding to each site name.
- species_vec: vector of length 134 corresponding to each species name.
- y: raw occurence data. 4 dimensional array of dimensions 52 (number of sites) * 7 (number of repeated surveys) * 17 (number of seasons) * 134 (number of species).

4. Missing data codes: NA

5. Specialized formats or other abbreviations used: NA

