#EXAMPLE ONE

#Let's explore the data visually 
plot(response~temperature, data=dat)

#add the line of best fit
fit1 = lm ("response~temperature", data=dat)
abline(fit1,col="red")

#view the spatial pattern of response, covariate and residuals
library(ggplot2)

x.range=range(x)
y.range=range(y)
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.01), y = seq(from = y.range[1], to = y.range[2], by = 0.01))  # expand points to grid
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE
xy=SpatialPoints(coords=cbind(x,y))
idw <- idw(formula = dat$temperature ~ 1, locations = xy, newdata = grd)  # apply idw model for the data
## [inverse distance weighted interpolation]
idw.output = as.data.frame(idw)  # output is defined as a data table
names(idw.output)[1:3] <- c("long", "lat", "var1.pred")  # give names to the modelled variables


ggplot() + geom_tile(data = idw.output, aes(x = long, y = lat, fill = var1.pred))+scale_fill_gradient2(low="blue", high="red", guide = guide_legend(title="temp")) +  geom_point(data = dat, aes(x = x, y = y, colour=as.factor(PresAbsresponse))) 
ggsave("Exploration_plot_example1.pdf")

ggplot() + geom_tile(data = idw.output, aes(x = long, y = lat, fill = var1.pred))+scale_fill_gradient2(low="blue", high="red", guide = guide_legend(title="temp")) 
ggsave("Temp_Map_example1.pdf")
