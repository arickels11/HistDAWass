library(HistDAWass)
library(waddR)
library(plyr)
library(dplyr)
library(purrr)

gas = `4_daily`
df = data.frame(gas)
summary(df)

# 48 after 5 removed, Ozone, Carbon monoxide, Sulfur dioxide, Nitrogen dioxide (NO2)
df_s1 = df[, c("State.Name","Parameter.Name", "Arithmetic.Mean")] %>%
  filter(!State.Name %in% c("Country Of Mexico","Alaska","Nebraska","West Virginia","Puerto Rico")) %>%
  filter(Parameter.Name == "Ozone")
df_s1c = df[, c("State.Name","Parameter.Name", "Arithmetic.Mean")] %>%
  filter(!State.Name %in% c("Country Of Mexico","Alaska","Nebraska","West Virginia","Puerto Rico")) %>%
  filter(Parameter.Name == "Carbon monoxide")
df_s1s = df[, c("State.Name","Parameter.Name", "Arithmetic.Mean")] %>%
  filter(!State.Name %in% c("Country Of Mexico","Alaska","Nebraska","West Virginia","Puerto Rico")) %>%
  #filter(Parameter.Name == "Sulfur dioxide")
  filter(Parameter.Name == "Sulfur dioxide", Arithmetic.Mean > 20) # outlier handling

df_s1n = df[, c("State.Name","Parameter.Name", "Arithmetic.Mean")] %>%
  filter(!State.Name %in% c("Country Of Mexico","Alaska","Nebraska","West Virginia","Puerto Rico")) %>%
  filter(Parameter.Name == "Nitrogen dioxide (NO2)")

# #SCALING
# df_s1["Arithmetic.Mean"] = scale(df_s1["Arithmetic.Mean"])
# df_s1c["Arithmetic.Mean"] = scale(df_s1c["Arithmetic.Mean"])
# df_s1s["Arithmetic.Mean"] = scale(df_s1s["Arithmetic.Mean"])
# df_s1n["Arithmetic.Mean"] = scale(df_s1n["Arithmetic.Mean"])


# group df by states
df_s2 = df_s1 %>%
  group_by(State.Name, Parameter.Name) %>% 
  group_map(~.x) %>%
  setNames(unique(sort(df_s1$State.Name)))
  
df_s2c = df_s1c %>%
  group_by(State.Name, Parameter.Name) %>% 
  group_map(~.x) %>%
  setNames(unique(sort(df_s1$State.Name)))
df_s2s = df_s1s %>%
  group_by(State.Name, Parameter.Name) %>% 
  group_map(~.x) %>%
  setNames(unique(sort(df_s1$State.Name)))
df_s2n = df_s1n %>%
  group_by(State.Name, Parameter.Name) %>% 
  group_map(~.x) %>%
  setNames(unique(sort(df_s1$State.Name)))


# transform to histogram distributions 
s2_d = map(df_s2, ~as.double(unlist(.x)))
o = lapply(s2_d, data2hist, type='regular', algo = "FixedQuantiles", qua = 30)

s2_d = map(df_s2c, ~as.double(unlist(.x)))
c = lapply(s2_d, data2hist, type='regular', algo = "FixedQuantiles", qua = 20)

s2_d = map(df_s2s, ~as.double(unlist(.x)))
s = lapply(s2_d, data2hist, type='regular', algo = "FixedQuantiles", qua = 8)

s2_d = map(df_s2n, ~as.double(unlist(.x)))
n = lapply(s2_d, data2hist, type='regular', algo = "FixedQuantiles", qua = 20)


hists = c(o,c,s,n)
states = unique(sort(df_s1$State.Name))

# Create HistDAWass matrix
m <- new("MatH", nrows = 48, ncols = 4, ListOfDist = hists,
         names.rows=states,
         names.cols=c("Ozone","Carbon_Monox","Sulf_Dioxide","Nit_Dioxide"),
         by.row = FALSE)



get.MatH.stats(m, stat='std')
cen_m <- Center.cell.MatH(m) # sets mean equal to 0

#Matrix plots
#plot(m, type = "HISTO", border = "blue") # plots a matrix of histograms
plot(m[10:20,c(1,4)], type = "DENS") # plots a matrix of densities
#plot(m, type = "BOXPLOT")


#breakout matrices by pollutant 

oz = (m[, "Ozone"])
carb = (m[, "Carbon_Monox"])
sulf = (m[, "Sulf_Dioxide"])
nit = (m[, "Nit_Dioxide"])


# Correlation, Covariance, Mean Histogram, Sum Histogram

WH.correlation(m)
covmat = WH.var.covar(m)

# vector of single pollutant
mean_hist = WH.vec.mean(nit)
plot(mean_hist)
get.MatH.stats(nit, stat = "max")
sum_hist = WH.vec.sum(oz)
plot(sum_hist)

# Matrix with 11 rows only
hists = c(nit[10:20,1],mean_hist)
states = unique(sort(df_s1$State.Name))

m <- new("MatH", nrows = 11, ncols = 4, ListOfDist = hists,
         names.rows=states,
         names.cols=c("Ozone","Carbon_Monox","Sulf_Dioxide","Nit_Dioxide"),
         by.row = FALSE)


WH.correlation(cen_m)
show(m)

me = WH.vec.mean(nit)
plot(me,type="DENS")
su = WH.vec.sum(oz)
plot(su)

plot(oz, type="HISTO")
plot(m["Florida", "Ozone"], type="HISTO")

#distances between states 
WassSqDistH(object1 = m["Hawaii", "Nit_Dioxide"]@M[[1]], object2 = m["Georgia", "Nit_Dioxide"]@M[[1]], details = TRUE)
WassSqDistH(object1 = m["Illinois", "Nit_Dioxide"]@M[[1]], object2 = m["Georgia", "Nit_Dioxide"]@M[[1]], details = TRUE)
WassSqDistH(object1 = m["Kansas", "Carbon_Monox"]@M[[1]], object2 = m["Maine", "Carbon_Monox"]@M[[1]], details = TRUE)
WassSqDistH(object1 = m["Georgia", "Carbon_Monox"]@M[[1]], object2 = m["Maine", "Carbon_Monox"]@M[[1]], details = TRUE)



# k means clustering
km = WH_adaptive.kmeans(x = oz, k = 5, rep = 10,
                              simplify = TRUE, qua =60, standardize = TRUE)

km$quality
plot(km$proto, type='HISTO')

#km$IDX
km$proto
table(km$IDX)
#oz .84 with 5
#carb .81 with 5, .877 with 8
#sulf .96 with 5 clusters
#nit .85 with 5 clusters

#all 4: .85 with 8, .77 with 5

#print out cluster assignment of each state
#clusters = bk$IDX
clusters = km$IDX

for (i in 1:length(clusters)) {
  (print(states[[i]]))
  (print(clusters[[i]]))
}

#register makes all p values the same (cdf's)
reg = registerMH(m)
reg

plot(reg[1:10], type = "DENS", border = "blue")

# fuzzy cmeans (gives weight to each cluster)

fc =  WH_adaptive_fcmeans(
  x = oz, k = 5, schema=4, m = 1.5, rep = 3, simplify = TRUE,
  qua = 10, standardize = TRUE, init.weights = "EQUAL", weight.sys = "PROD"
)

plot(fc$solution$proto[1], type="DENS")

fc$solution$IDX
fc$solution$membership
fc$quality
fc$solution$proto


# agglom clustering

hc = WH_hclust(
  oz,
  simplify = FALSE,
  qua = 10,
  standardize = FALSE,
  distance = "WDIST",
  method = "complete"
)

plot(hc) # it plots the dendrogram
cut = cutree(hc, k = 6) # it returns the labels for 5 clusters
View(cut)
cut
hc

# non Wass clustering

library(stats)
hce = hclust(dist(x=s2_d, method="euclidean"))
plot(hce)
memb = cutree(hce, k = 5)
memb

hc1 = stats.hclust(dist(cent)^2, method = "cen", members = table(memb))
opar = par(mfrow = c(1, 2))
plot(hce,  labels = FALSE, hang = -1, main = "Original Tree")
plot(hc1, labels = FALSE, hang = -1, main = "Re-start from 10 clusters")
par(opar)





# Correlation between variables (not state-wise)

df_g = filter(df, !State.Name %in% c("Country Of Mexico","Alaska","Nebraska","West Virginia","Puerto Rico"))
#unique(df_g$State.Name)
df_gas = df_g[, c("Parameter.Name", "Arithmetic.Mean")] 

df_gas_list = df_gas %>%
  filter(!is.na("Arithmetic.Mean")) %>%
  group_by(Parameter.Name) %>% 
  group_map(~.x) %>%
  setNames(unique(sort(df_gas$Parameter.Name)))

  
gas_list = map(df_gas_list, ~as.double(unlist(.x)))
gas_dist = lapply(gas_list, data2hist, type='regular', algo = "FixedQuantiles", qua = 20) 

#gas_dist = lapply(gas_list, scale) # scale data
#gas_dist = lapply(gas_dist, data2hist) # scale data




#W-squared distances between pollutant distributions
WassSqDistH(object1 = gas_dist[["Ozone"]], object2 = gas_dist[["Carbon monoxide"]], details = TRUE)
WassSqDistH(object1 = gas_dist[["Ozone"]], object2 = gas_dist[["Nitrogen dioxide (NO2)"]], details = TRUE)
WassSqDistH(object1 = gas_dist[["Ozone"]], object2 = gas_dist[["Sulfur dioxide"]], details = TRUE)
WassSqDistH(object1 = gas_dist[["Nitrogen dioxide (NO2)"]], object2 = gas_dist[["Carbon monoxide"]], details = TRUE)
WassSqDistH(object1 = gas_dist[["Nitrogen dioxide (NO2)"]], object2 = gas_dist[["Sulfur dioxide"]], details = TRUE)
WassSqDistH(object1 = gas_dist[["Carbon monoxide"]], object2 = gas_dist[["Sulfur dioxide"]], details = TRUE)

plot(gas_dist[["Ozone"]], type = "CDF", col = "green", border = "black")
plot(gas_dist[["Ozone"]], type = "QF", col = "green", border = "black")


plot(gas_dist[["Nitrogen dioxide (NO2)"]], type = "DENS", col = "blue", border = "black")
mean_hist = WH.vec.mean(nit)
plot(mean_hist, type="HISTO")
plot(m[,4],type = "HISTO")

#rQQ is correlation coefficient between quantile functions

m_all <- new("MatH", nrows = 1, ncols = 4, ListOfDist = gas_dist,
             names.cols=c("Carbon_Monox","Nit_Dioxide","Ozone","Sulf_Dioxide"),
             by.row = FALSE)


get.MatH.stats(m_all)
plot(m_all, type="DENS")
WH.correlation(m_all)

#Matrix plots
plot(m_all[["Ozone"]], type = "DENS", border = "blue") # plots a matrix of densities

plot(gas_dist[["Ozone"]], type = "DENS", col = "yellow", border = "black")
plot(gas_dist[["Carbon monoxide"]], type = "DENS", col = "blue", border = "black")
plot(gas_dist[["Nitrogen dioxide (NO2)"]], type = "DENS", col = "red", border = "black")
plot(gas_dist[["Sulfur dioxide"]], type = "DENS", col = "green", border = "black")

#updated for outliers

# sulfur dioxide: 425 / 151K over 20
# carbon monox: 4 / 86K over 2.2 (all puerto rico)

sulf = gas_list$`Sulfur dioxide`[gas_list$`Sulfur dioxide` < 10]
sulf_dist = data2hist(sulf)
car = gas_list$`Carbon monoxide`[gas_list$`Carbon monoxide` < 2.2]
car_dist = data2hist(car)

oz_dist = gas_dist[["Ozone"]]
nit_dist = gas_dist[["Nitrogen dioxide (NO2)"]]

hists_new = c(oz_dist, nit_dist, sulf_dist, car_dist)

m_all_new <- new("MatH", nrows = 1, ncols = 4, ListOfDist = hists_new,
             names.cols=c("Ozone","Nit_Dioxide","Carbon_Monox","Sulf_Dioxide"),
             by.row = FALSE)


get.MatH.stats(m_all_new)
plot(m_all_new, type="DENS")

#new W-squared distances - all above .84
WassSqDistH(object1 = oz_dist, object2 = car_dist, details = TRUE)
WassSqDistH(object1 = oz_dist, object2 = nit_dist, details = TRUE)
WassSqDistH(object1 = oz_dist, object2 = sulf_dist, details = TRUE)
WassSqDistH(object1 = nit_dist, object2 = car_dist, details = TRUE)
WassSqDistH(object1 = nit_dist, object2 = sulf_dist, details = TRUE)
WassSqDistH(object1 = car_dist, object2 = sulf_dist, details = TRUE)

WH.correlation(m_all_new)
WH.var.covar(m_all_new)

plot(oz_dist, type="QF")
plot(car_dist, type="QF")
plot(nit_dist, type="QF")
plot(sulf_dist, type="QF")



# Cities in California
df_city = df[, c("State.Name", "City.Name", "Parameter.Name", "Arithmetic.Mean")] %>%
  filter(State.Name %in% c("California")) %>%
  filter(Parameter.Name == "Nitrogen dioxide (NO2)")

#df_s1n["Arithmetic.Mean"] = scale(df_s1n["Arithmetic.Mean"])

df_cities = df_city %>%
  group_by(State.Name, City.Name, Parameter.Name) %>% 
  group_map(~.x) %>%
  setNames(unique(sort(df_s1$City.Name)))

cities = map(df_cities, ~as.double(unlist(.x)))
cities_hist = lapply(cities, data2hist)


cities = unique(sort(df_city$City.Name))

m_cities <- new("MatH", nrows = 77, ncols = 1, ListOfDist = hists,
         names.rows=cities,
         names.cols=c("Nit Diox"),
         by.row = FALSE)

plot(m_cities, type='DENS')

# agglom clustering of cities

hc = WH_hclust(
  m_cities,
  simplify = FALSE,
  qua = 10,
  standardize = FALSE,
  distance = "WDIST",
  method = "complete"
)

plot(hc) # it plots the dendrogram
cut = cutree(hc, k = 5) # it returns the labels for 5 clusters
cut

dend = as.dendrogram(cut)
plot(dend)

colors(dend) = 2:4
labels_colors(dend)
plot(dend)
