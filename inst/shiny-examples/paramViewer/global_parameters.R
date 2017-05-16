## Names of all the allowable strains. Can change this to actual strain names
## if preferred
exposure_strains <- c("A","B","C","D","E",
                      "F","G","H","I","J")
weak_types <- c("Infection 1"="infection1","Infection 2"="infection2", "Vaccine 1"="vacc1",
                "Vaccine 2"="vacc2","Adjuvanted 1"="adj1","Adjuvanted 2"="adj2")
strong_types <- c("Infection"="infection","Vaccine"="vacc","Adjuvanted"="adj")

max_mu <- 20
max_tp <- 100
max_ts <- 100
min_m <- -10
