# Spatial Mapping; part of parameter evaluation
# Spatial data sets come from two sources:
# (1) Tiger/Line & US Census Bureau -- shapefiles
# (2) Google Maps (is the best, true that, double true)

rm(list=ls())


require(maptools)
require(gpclib) #dependency for maptools, wasn't automatically loaded the first go-around, so I'm adding this line just in case
require(stringr)
library(RColorBrewer)
require(ggplot2)
require(scales)
gpclibPermit()

get.centroids<-function(batch.data, directory){
  # batch.data is the name of a file that contains:
  #       the name of the county, if applicable, written as stored in the shape file
  #       the name of the state, written as in the shape file name
  shapes<-list.files(directory, pattern=".shp$", full.names=TRUE, recursive=TRUE, include.dirs=TRUE) #search the directory for the shape files
  #print(shapes)
  lat=c()
  long=c()
  for(i in 1:length(batch.data$State)){
    if(is.na(batch.data$County[i])==1){ #if there is no county name given, get the centroid for the state as a whole
      
      s.name<-batch.data$State[i] #what state are we looking at?
      regex<-paste('state10_', s.name, sep="") #build the fact that we're looking for a state, and the name of the state, into a regex
      s<-str_detect(shapes, pattern=regex) #which index matches our pattern? (returns a vector of T/F)
      j=which(s == TRUE) #pull out the index that's true... if we've done it right, there will only be one
      state<-readShapePoly(shapes[j]) #read the file from the name we found above
      coord=coordinates(state) #get centroid
      long<-c(long, coord[1])
      lat<-c(lat, coord[2])}
      
    else{
      s.name<-batch.data$State[i] #what state are we looking at?
      regex<-paste('county10_', s.name, sep="") #build the fact that we're looking for a state, and the name of the state, into a regex
      s<-str_detect(shapes, pattern=regex) #which index matches our pattern? (returns a vector of T/F)
      j=which(s == TRUE) #pull out the index that's true... if we've done it right, there will only be one
      print(j)
      print(shapes[j])
      state.counties<-readShapePoly(shapes[j]) #read the file from the name we found above
      c.name<-as.character(batch.data$County[i]) #can't pass the name as a factor; conver to character
      final.shape<-state.counties[state.counties$NAME10==c.name,] #get the data for county of interest only
      coord=coordinates(final.shape) #get centroid
      long<-c(long, coord[1])
      lat<-c(lat, coord[2])}
    }
  coords<-cbind(batch.data, long, lat)
  return(coords)
}


# get data
bd<-read.csv("/Users/kelly/Dropbox/2014&Older/ModelFitting/LitSearches_Data/Counties_States.csv")
dir=c("/Users/kelly/Dropbox/2014&Older/ModelFitting/ShapeFiles")

bd2<-read.csv("/Users/kelly/Dropbox/2014&Older/ModelFitting/LitSearches_Data/Counties_States_2.csv")

dir2="/Users/kelly/Dropbox/2014&Older/ModelFitting/ShapeFiles/tl_2010_13_county10_Georgia/"

#output data
new.data<-get.centroids(bd, dir)
write.csv(new.data, file='relevant_centroids.csv')

new.data2<-get.centroids(bd2, dir)
write.csv(new.data2, file='relevant_centroids_2.csv')
# these data were reogranized in Excel to make files of point coords for the different data types (see below)

# get shapefiles
USA<-readShapePoly("/Users/kelly/Dropbox/2014&Older/ModelFitting/ShapeFiles/tl_2010_us_state10/tl_2010_us_state10.shp")

# make an inset map to give context to Southeast US
#lwr.48 is intentially missing Alabama
lwr.48<-c('Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Florida','Georgia','Idaho','Illinois',
          'Indiana','Iowa','Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts','Michigan','Minnesota','Mississippi',
          'Missouri','Montana','Nebraska','Nevada','New Hampshire','New Jersey','New Mexico','New York','North Carolina','North Dakota','Ohio',
          'Oklahoma','Oregon','Pennsylvania','Rhode Island','South Carolina','South Dakota','Tennessee',
          'Texas','Utah','Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming')
lwr.48.shp<-USA[USA$NAME10=='Alabama',] #start dataframe w/Alabama
for(state in lwr.48){
  shape<-USA[USA$NAME10==state,]
  lwr.48.shp<-rbind(lwr.48.shp, shape)
}
par(mfrow=c(1,1))
plot(lwr.48.shp)

## to make a plot of just the points in the southeast
south.east<-c('Florida', 'North Carolina', 'South Carolina')
se<-USA[USA$NAME10=='Georgia',]
for(state in south.east){
  shape<-USA[USA$NAME10==state,]
  se<-rbind(se, shape)
}

# get outlines for surrounding states to give map some more context
se.plus<-se # give se a new name
south.east.extra<-c('Virginia', 'Kentucky', 'Tennessee', 'Alabama')
for(state in south.east.extra){
  shape<-USA[USA$NAME10==state,]
  se.plus<-rbind(se.plus, shape)
}

# thin the shape files for faster/easier plotting
thin.48<-thinnedSpatialPoly(lwr.48.shp, tolerance=0.05)
thin.se<-thinnedSpatialPoly(se, tolerance=0.05)
thin.se.plus<-thinnedSpatialPoly(se.plus, tolerance=0.05)

# plot the lower 48 states and shade the southeast ones for which we have data in gray
pdf('~/Dropbox/2014&Older/ModelFitting/full_us_map.pdf', height=20/2.54, width=20/2.54)#, units='cm', res=300)
plot(thin.48, lty=3)
plot(thin.se, col='gray', add=TRUE)
dev.off()

# read in reorganized data points for each study type
# we won't actually use raccoon data... 
#racc.points<-read.csv('/Users/kelly/Dropbox/2014&Older/ModelFitting/LitSearches_Data/Racc_Coords.csv')
deer.points<-read.csv('/Users/kelly/Dropbox/2014&Older/ModelFitting/LitSearches_Data/Deer_Coords.csv')
tick.points<-read.csv('/Users/kelly/Dropbox/2014&Older/ModelFitting/LitSearches_Data/Tick_Coords.csv')

# the deer point that is marked PCR/IFA is an infection prevalence result
deer.ab.points<-deer.points[deer.points$Technique=='IFA',]
deer.pcr.points<-deer.points[deer.points$Technique!='IFA',]

# lots of different color options...
#col.pal=c(rev(brewer.pal(3,'Greens'))[1], rev(brewer.pal(3,'Blues'))[1], rev(brewer.pal(3,'Reds'))[1])
#print.col=c('#ef8a62', '#67a9cf', '#ffffbf')
print.col=c('violetred3', 'blue1', 'yellow3')

# lots of different plot formatting options
#pdf('southeast_map_for_presentation_whiteBG.pdf', width=10, height=7)
#png('southeast_map_for_presentation_whiteBG.png', width=20, height=14, unit='cm', res=300)
#pdf('southeast_map_for_presentation_whiteBG_extra_states.pdf', width=20/2.54, height=20/2.54)#, unit='cm', res=300)
png('southeast_map_for_presentation_whiteBG_extra_states.png', width=20, height=20, unit='cm', res=300)

#par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="white", mar=c(0.5,0.5,0.5,0.5))
#plot(se, lwd=3)

par(mfrow=c(1,1), mar=c(2,2,2,4))
plot(thin.se.plus, lty=3) # plot the extended states with dotted borders
plot(thin.se, lwd=3, add=T) # add the SE states from which we have data with thick solid borders

# add the points (many different styles listed below)

# colors from presentation
#points(jitter(tick.points$Longitude, amount=0), jitter(tick.points$Latitude, amount=0), pch=16, col=alpha(rev(brewer.pal(3,'Greens'))[1], 0.7), cex=2)
#points(jitter(deer.points$Longitude, amount=0), jitter(deer.points$Latitude, amount=0), pch=16, col=alpha(rev(brewer.pal(3,'Blues'))[1], 0.7), cex=2)

# print and colorblind safe colors
points(jitter(tick.points$Longitude, amount=0), jitter(tick.points$Latitude, amount=0), pch=21, bg=alpha(print.col[1], 0.5), cex=2)
#points(jitter(deer.points$Longitude, amount=0), jitter(deer.points$Latitude, amount=0), pch=16, col=alpha(print.col[2], 0.5), cex=2)

# can split deer up by PCR vs. AB, but there's not much difference in sites -- most studies did both.
points(jitter(deer.ab.points$Longitude, amount=0), jitter(deer.ab.points$Latitude, amount=0), pch=21, bg=alpha(print.col[2], 0.5), cex=2)
points(jitter(deer.pcr.points$Longitude, amount=0), jitter(deer.pcr.points$Latitude, amount=0), pch=21, bg=alpha(print.col[3], 0.5), cex=2)

#points(jitter(racc.points$Longitude, amount=0), jitter(racc.points$Latitude, amount=0), pch=16, col=alpha(rev(brewer.pal(3,'Reds'))[1], 0.5), cex=2)
#use this legend to plot 4 states only and properly position legend
#legend(x=-80, y=32, legend=c('Adult Ticks', 'Deer'), col=c(alpha(col.pal[1:2], 1)), bty='n', pch=16, cex=2)
#use this legend to plot the extra se states and properly position legend
legend(x=-80, y=33, legend=c('Adult Ticks', 'Deer', 'Deer AB'), pch=21, col='black', pt.bg=c(alpha(print.col, 1)), bty='n', cex=1.5)
dev.off()

## NOTE: the complete figure was created by merging the lower 48 states map and the data point map
## in powerpoint (inkscape would be an option with PDF format)

## Various maps of Texas with field site locations marked
setwd('~/Dropbox/')
texas<-USA[USA$NAME10=='Texas',]
png('texas_outline_CLNWR_white_bg_small.png', width=10, height=10, unit='cm', res=300)
plot(texas, lwd=3)
points(-94.1, 32.6, pch=16, col='blue', cex=3)
dev.off()

setwd('~/Dropbox/ModelFitting/FutureProof/Figures for Presentations/')
texas<-USA[USA$NAME10=='Texas',]
png('texas_outline_black_bg_small.png', width=10, height=10, unit='cm', res=300)
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="red", mar=c(0.5,0.5,0.5,0.5))
plot(texas, lwd=3)
dev.off()

png('texas_outline_GEWMA_black_bg_small.png', width=10, height=10, unit='cm', res=300)
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="red", mar=c(0.5,0.5,0.5,0.5))
plot(texas, lwd=3)
points(-95.9, 31.9, pch=16, col='blue', cex=3)
dev.off()

png('texas_outline_GEWMA_white_bg_small.png', width=10, height=10, unit='cm', res=300)
par(mar=c(0.5,0.5,0.5,0.5))
plot(texas, lwd=3)
points(-95.9, 31.9, pch=16, col='blue', cex=3)
dev.off()

png('texas_outline_CLNWR_black_bg_small.png', width=10, height=10, unit='cm', res=300)
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="red", mar=c(0.5,0.5,0.5,0.5))
plot(texas, lwd=3)
points(-94.1, 32.6, pch=16, col='blue', cex=3)
dev.off()

png('texas_outline_Austin_black_bg_small.png', width=10, height=10, unit='cm', res=300)
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="red", mar=c(0.5,0.5,0.5,0.5))
plot(texas, lwd=3)
points(-97.7, 30.3, pch=16, col='blue', cex=3)
dev.off()

# -------------------------------------------------------------------------------------
# -- Various notes and test scripts on working with shapefiles and color interpolation

# Unclear why I have the eastern states, but here they are:
eastern.states<-c('Vermont', 'New Hampshire', 'Rhode Island', 'Massachusetts', 'New York', 'Pennsylvania', 'Ohio', 'Texas', 'Oklahoma', 'Kansas', 'Nebraska', 'South Dakota', 'North Dakota', 'Michigan', 'Indiana', 'Illinois', 'Wisconsin', 'Minnesota', 'Iowa', 'Missouri', 'Arkansas', 'Mississippi', 'Alabama', 'Georgia', 'Florida', 'Tennessee', 'Kentucky', 'Virginia', 'West Virginia', 'North Carolina', 'South Carolina', 'Maryland', 'New Jersey', 'Delaware', 'Connecticut', 'District of Columbia', 'Louisiana')
#Maine is missing from the above vector because we use it to start off the data frame
east<-USA[USA$NAME10=='Maine',]
for(state in eastern.states){
  shape<-USA[USA$NAME10==state,]
  east<-rbind(east, shape)}

plot(east) #makes a nice pretty map of the eastern US!
# ggplot(east) #don't try this; it hangs

adult.prev<-read.csv('/Users/kellypierce/Dropbox/ModelFitting/Adult_Tick_Inf_Prev.csv')
deer.prev<-read.csv('/Users/kellypierce/Dropbox/ModelFitting/Deer_Inf_Prev.csv')
# all lat/longs should be filled in now...
#ap<-na.omit(adult.prev) #remove ones for which there is not yet a lat/long
#dp<-na.omit(deer.prev)

prev.cols<-colorRamp(colors=c('#BFD2E6'), space=c('rgb')) #get an interpolation function between these two colors -- we can call prev.cols as a function because it is!
cs<-prev.cols(adult.prev$Prevalence)

breaks=seq(from=0, to=1, by=0.1) #find breaking points between 0 and 100% prevalence in 10% increments
bins<-cut(adult.prev$Prevalence, breaks=breaks, include.lowest=TRUE) #returns the bin ID for each value in adult.prev$Prevalence
adult.prev<-cbind(adult.prev, bins) #add the bin values into the data set

# finally! a way to get a color palette with more colors in the same hue: http://www.perbang.dk/rgbgradient/
# this gives you opaque colors. want them transparent?
# Transparency is scaled 0-255 (nothing to fully opaque). To get an intermediate transparency, convert 0-225 to a 0-100 percentage (e.g., 50% = 127), and convert that number to a HEX value.
# Hex conversions can be made here: http://www.binaryhexconverter.com/decimal-to-hex-converter
# To inform R that you want transparent colors, append the two-digit HEX value to the beginning of your HEX color code, but after the '#':
#       '#7F------' will give you 50% transparency in whatever color you choose (127 in HEX is 4B)
# this one is ordered from dark to light...
my.blues.pal<-c('#7F0124A2','#7F1131A7','#7F213FAC','#7F324DB1','#7F425BB6','#7F5269BC','#7F6377C1','#7F7385C6','#7F8393CB','#7F94A1D1')
# flip it.
my.blues.pal<-rev(my.blues.pal)
# this one is already ordered light to dark
my.reds.pal<-c('#7FCD9493','#7FC88482','#7FC37472','#7FBE6462','#7FB95452','#7FB54441','#7FB03431','#7FAB2421','#7FA61411','#7FA20401')

plot(east)
points(jitter(adult.prev$Longitude), jitter(adult.prev$Latitude), pch=16, col=my.blues.pal[bins])
points(jitter(deer.prev$Longitude), jitter(deer.prev$Latitude), pch=16, col=my.reds.pal[bins])

# one example of how to get the colors on the plog... rgb() allows you to select the intensity and transparency of either red, green, or blue points
#points(jitter(deer.prev$Longitude), jitter(deer.prev$Latitude), pch=16, col=rgb(0,0,100,50,maxColorValue=225))

# another example of getting the colors... this works nicely (bins defined as above), but I don't know how to get colors other than green!
#colors<-terrain.colors(n=length(breaks), alpha=0.5)[bins] # get a color for each bin, and make a vector that contains the colors that belong to each bin in bins

#generic syntax for plotting a county that is part of a state's shapefile
dk<-GA.co[GA.co$NAME10=='DeKalb',]
plot(dk)
fl<-GA.co[GA.co$NAME10=='Fulton',]
plot(fl)

atl<-rbind(dk, fl) #this merges the shapes, kind of
plot(atl) #see what I mean?
atl.coords<-coordinates(atl) #you get two centroids
points(atl.coords)

#get the polygon centroid; from the sp package required by maptools
coords<-coordinates(GA)

#this is all just in a file now...
#in.cos<-c('Crawford', 'Orange', 'Pike', 'Perry', 'Spencer', 'Dubois', 'Warrick')
#ga.cos<-c('Ogelthorpe')
#ia.cos<-c('Wapello')
#ok.cos<-c('Payne')
#oh.cos<-c('Adams', 'Gallia', 'Jackson', 'Lawrence')
#ms.cos<-c('Marion')
#nc.cos<-c('Montgomery')
#va.cos<-c('Roanoke')
#wv.cos<-c('Wood')
# la.cos<-c('Caldwell', 'Franklin') these need to be joined together and a joint centroid taken 

'''

data.GA<-read.csv("/Users/kellypierce/Dropbox/ModelFitting/Partial_Adult_Prev_Geo.csv")

I already wrote it all out, so why delete it?

# State ShapeFiles -- these load properly and are the correct states!
AL<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_01_state10_Alabama/tl_2010_01_state10.shp")
FL<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_12_state10_Florida/tl_2010_12_state10.shp")
GA<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_13_state10_Georgia/tl_2010_13_state10.shp")
IL<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_17_state10_Illinois/tl_2010_17_state10.shp")
KY<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_21_state10_Kentucky/tl_2010_21_state10.shp")
MD<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_24_state10_Maryland/tl_2010_24_state10.shp")
MA<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_25_state10_Massachusetts/tl_2010_25_state10.shp")
MO<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_29_state10_Missouri/tl_2010_29_state10.shp")
NC<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_37_state10_NorthCarolina/tl_2010_37_state10.shp")
TX<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_48_state10_Texas/tl_2010_48_state10.shp")

# Counties ShapeFiles -- these load properly and are the correct states with counties!
FL.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_12_county10_Florida/tl_2010_12_county10.shp")
GA.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_13_county10_Georgia/tl_2010_13_county10.shp")
IN.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_18_county10_Indiana/tl_2010_18_county10.shp")
IA.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_19_county10_Iowa/tl_2010_19_county10.shp")
LA.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_22_county10_Louisiana/tl_2010_22_county10.shp")
MS.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_28_county10_Mississippi/tl_2010_28_county10.shp")
NC.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_37_county10_NorthCarolina/tl_2010_37_county10.shp")
OH.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_39_county10_Ohio/tl_2010_39_county10.shp")
OK.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_40_county10_Oklahoma/tl_2010_40_county10.shp")
VA.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_51_county10_Virginia/tl_2010_51_county10.shp")
WV.co<-readShapePoly("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/tl_2010_54_county10_WestVirginia/tl_2010_54_county10.shp")
'''

'''
snippets of test code:

#this line finds all the files
shapes<-list.files("/Users/kellypierce/Dropbox/ModelFitting/ShapeFiles/", pattern=".shp$", full.names=TRUE, recursive=TRUE, include.dirs=TRUE)

#generic syntax for getting all the file names, finding which contain the desired pattern, and getting the path for that file
texas<-str_detect(shapes, pattern='Texas')
i=which(texas == TRUE)
name=shapes[i]
test<-readShapePoly(shapes[1])

x11()
plot(GA) #plot the outline of the shape
points(coords, pch=16)
points(as.numeric(data.GA$Long), as.numeric(data.GA$Lat), pch=16, col='red')


'''
