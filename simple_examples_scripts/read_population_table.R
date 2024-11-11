# Read population table

library(relsurv)

pop.data.male <-
  read.table("norwegian_population_tables/NOR.mltper_1x1.txt",
             skip = 2,
             header = T) # Life table for male in Norway
pop.data.female <-
  read.table("norwegian_population_tables/NOR.fltper_1x1.txt",
             skip = 2,
             header = T) # Life table for female in Norway
nortab <-
  transrate.hmd(male = "norwegian_population_tables/NOR.mltper_1x1.txt", female = "norwegian_population_tables/NOR.fltper_1x1.txt")
