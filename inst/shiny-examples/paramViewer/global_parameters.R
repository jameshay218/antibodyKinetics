## Names of all the allowable strains. Can change this to actual strain names
## if preferred
exposure_strains <- c("A","B","C","D","E",
                      "F","G","H","I","J")
pairs <- 
weak_types <- c("Infection 1"="infection1","Infection 2"="infection2", "Vaccine 1"="vacc1",
                "Vaccine 2"="vacc2","Adjuvanted 1"="adj1","Adjuvanted 2"="adj2")
strong_types <- c("Infection"="infection","Vaccine"="vacc","Adjuvanted"="adj")

max_mu <- 8
min_mu <- 2

max_dp <- 0.7
min_dp <- 0.4

max_tp <- 17
min_tp <- 10

max_ts <- 15
min_ts <- 5

max_m <- -3
min_m <- -8
