## Preamble ####
# The following script is a "port" of Borgatti, Everett, & Johnson (2nd
#  edition) to R. (Well, at least it is a port of most of it.)
# In a couple locations, I have used fread() to open ".dat" files, which is
#  a common network file format. I have also included a couple functions to
#  read basic UCINET files, which are also a common format. However, I think
#  that all the data are also included as .CSV files; I ported the data from
#  various sources.
#
# This isn't in "package" form because I want you to have the code, and
#  not just the output. (I don't know if you knew this, but you can
#  get the code for any function by typing the function name without the
#  parentheses at the console, and then pressing "Enter". But I wanted
#  to make it easier for you.)
#
# This was prepared/updated for use by GDAT 622 students at
#  St. John Fisher College; hence, the comments about what is
#  known or not known apply to that audience.
#

## Preparation ####

#
# There are several other potentially useful packages, including:
#   CINNA - 
#   gplot -
#   ggnetwork - ggplot-like plotting for networks
#   graphlayouts - Some useful plotting stuff for igraph
#   keyplayer - advanced centality methods
#   multiplex - for multplex graphs
#   ndtv - Network Dynamic Temporal Visualization
#   qgraph
#   RSiena - particularly useful for longitudinal networks; see Snijders' 
#            Siena home page
#   sand - from a more mathematical text
#   tnet - miscellaneous advanced analyses
#   vegan - Some useful network and clustering functions
# They are mostly beyond the scope of this course, however, and so are not
#  included here.

{
  c("alphahull",     # To calculate the convex hull
    "ca",            # Correspondence analysis
    "conflicted",    # To deal with conflicting function names
    # I've had some strangeness with this
    #  script. I suspect package:conflicted,
    #  but I don't yet know for sure.
    "data.table",    # Fast data input/output
    "dplyr",         # This is all of tidyverse that gets used here
    "dtplyr",        # dplyr syntax with a data.table backend
    "here",          # To find/store files w/o setwd() and getwd()
    "igraph",        # Basic network tools; we'll use statnet mostly
    "igraphdata",    # Some useful datasets
    "intergraph",    # Translate between igraph and statnet formats
    "lmPerm",        # To do permutation tests
    "statnet",      # A suite of network tools, including ERGM and more
    "xlsx"
  ) -> package_names
  
  for (package_name in package_names) {
    if (!is.element(package_name, installed.packages()[, 1])) {
      install.packages(package_name,
                       repos = "http://cran.mtu.edu/")
      # An alternate, just in case.
      #                      repos="http://lib.stat.cmu.edu/R/CRAN")
    }
    library(
      package_name,
      character.only = TRUE,
      quietly = TRUE,
      verbose = FALSE
    )
  }
  rm(list = c("package_name", "package_names"))
}

# Because I like these options:
options(show.signif.stars = FALSE)
options(digits = 4)

## Reading UCINET files ####
# package:rucinet used to read UCINET files, but UCINET has changed their
#  format, so that the rucinet package no longer works. However, I found some
#  other code, modified it - so we have access to the read.ucinet.header() 
#  and read.ucinet() functions - and include it here.
#
# Original code written by Christian Steglich,
#  https://sites.google.com/site/ucinetsoftware/document/faq/connectingwithr

read.ucinet.header <- function(filename) {
  # function for reading UCINET header files (recent versions only)
  # This only works for "single level" files. (In spreadsheet parlance, 
  # each workbook can have only one sheet.)
  UCINET.header <- file(paste(filename,".##h",sep=''),"rb")
  ignore <- readBin(UCINET.header,what="int",size=1)
  headerversion <- paste(
    rawToChar(readBin(UCINET.header,what="raw",size=1)),
    rawToChar(readBin(UCINET.header,what="raw",size=1)),
    rawToChar(readBin(UCINET.header,what="raw",size=1)),
    rawToChar(readBin(UCINET.header,what="raw",size=1)),
    rawToChar(readBin(UCINET.header,what="raw",size=1)),
    sep='')
  if (!(headerversion %in% c('DATE:','V6404'))) {
    close(UCINET.header)
    stop(paste('unknown header type; try more recent UCINET file types'))
  }
  year <- 2000+readBin(UCINET.header,what="int",size=2)
  month <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',
             'Sep','Oct','Nov','Dec')[readBin(UCINET.header,what="int",size=2)]
  day <- readBin(UCINET.header,what="int",size=2)
  dow <- c('Monday','Tuesday','Wednesday','Thursday','Friday',
           'Saturday','Sunday')[readBin(UCINET.header,what="int",size=2)]
  labtype <- readBin(UCINET.header,what="int",size=2)
  infile.dt <- c('nodt','bytedt','booleandt','shortintdt','worddt',
                 'smallintdt','longintdt',' singledt','realdt','doubledt',
                 'compdt','extendeddt','labeldt','setdt','stringdt','pointerdt',
                 'chardt','integerdt','nodelistdt','sparsedt','int64dt')[
                   readBin(UCINET.header,what="int",size=1)]
  ndim <- readBin(UCINET.header,what="int",size=2)
  if (headerversion=='V6404') {fct=2} else {fct=1}
  dims <- c(readBin(UCINET.header,what="int",size=2*fct),
            readBin(UCINET.header,what="int",size=2*fct))
  if (ndim==3) {
    dims[3] <- readBin(UCINET.header,what="int",size=2*fct)
  }
  if (!(ndim==2|ndim==3&dims[3]==1)) {
    close(UCINET.header)
    stop(paste('UCINET file with',dims[3],'levels; please convert separately'))
  }
  t.length <- readBin(UCINET.header,what="int",size=1)
  if (t.length>0){
    titl <- sapply(1:t.length, function(i){
      rawToChar(readBin(UCINET.header,what="raw",size=1))
    })
    titl <- paste(titl,collapse='')
  } else {titl <- ''}
  haslab <- c(readBin(UCINET.header,what="logical",size=1),
              readBin(UCINET.header,what="logical",size=1))
  if (ndim==3) {
    haslab[3] <- readBin(UCINET.header,what="logical",size=1)
  }
  dim.labels <- list()
  for (arr.dim in 1:length(dims)) {
    if (haslab[arr.dim]) {
      dim.labels[[arr.dim]] <- rep(NA,dims[arr.dim])
      for (i in 1:dims[arr.dim]) {
        lab <- ''
        lablen <- readBin(UCINET.header,what="int",size=2)
        for (let in 1:lablen) {
          lab <- paste(lab,
                       rawToChar(readBin(UCINET.header,what="raw",size=1)),
                       sep='')
        }
        dim.labels[[arr.dim]][i] <- lab
      }
    }}
  close(UCINET.header)
  if (ndim==3&dims[3]==1) {
    titl <- dim.labels[[3]][1]
    warning(paste('UCINET file with one level; level name "',
                  titl,'" treated as network name',sep=''))
    ndim <- 2
    dims <- dims[1:2]
    haslab <- haslab[1:2]
    dim.labels <- dim.labels[1:2]
  }
  return(list(
    headerversion=headerversion,
    date=paste(dow,paste(day,month,year,sep='-')),
    labtype=labtype,
    infile.dt=infile.dt,
    ndim=ndim,
    dims=dims,
    title=titl,
    haslab=haslab,
    dim.labels=dim.labels
  ))
}

read.ucinet <- function(filename) {
  # function for reading UCINET data files (recent versions only)
  # filename = UCINET filename (without ## extension)
  # begin of main function code:
  header <- read.ucinet.header(filename)
  UCINET.data <- file(paste(filename,".##d",sep=''),"rb")
  thedata <- c()
  for (i in 1:(header$dims[1]*header$dims[2]))
    thedata[i] <- readBin(UCINET.data,what="numeric",size=4,endian='little')
  close(UCINET.data)
  mat <- matrix(thedata,nr=header$dims[2],nc=header$dims[1],
                dimnames=header$dim.labels[c(2,1)],byrow=TRUE)
  # put additional info from header file on matrix
  if (header$title!='') {attr(mat,'title') <- header$title}
  attr(mat,'date') <- header$date
  #attr(mat,'labtype') <- header$labtype
  #attr(mat,'infile.dt') <- header$infile.dt
  return(mat)
}

## Chapter 1 ####
# No code

## Chapter 2 ####

# Figure 2.1
network.initialize(n = 5,
                   directed = TRUE) -> F2.1
LETTERS[1:5] -> network.vertex.names(F2.1)
network::add.edges(F2.1,
                   #          tail = c("A", "B", "C", "D", "D"),
                   tail = c(1, 2, 3, 4, 4),
                   head = c(2, 3, 4, 5, 1))
#          head = c("B", "C", "D", "A", "E"))

# A Note about plotting:
# Generic network plotting has a random aspect to it, so setting the seed
#  helps to make it reproducible. There are a number of ways to layout the
#  nodes on a network, and we'll look at some of those, but network plotting
#  is still a bit of an art and a bit of luck.
# Tcl/Tk can be used to set the position of nodes by hand, but OMG, what a
#   pain!
#
# All that is the long explanation for why the plots here will not necessarily
#  look exactly like those in the text.
# 
# More on plotting will happen in Chapter 7.

set.seed(42)
plot.network(F2.1,
             label = network.vertex.names(F2.1),
             arrowhead.cex = 2,
             edge.lwd = 2,
             vertex.col = 8,
             vertex.cex = 3,
             vertex.border = 0)

# The data, and descriptions of those data, for Figure 2.2 can be found at:
#   http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm#wiring
# Getting these data into R is not as much fun as one might like. Welcome to
#  the network version of data wrangling!

# Figure 2.2 common data for (a) and (b). This is a common thing to have to do
#  to get data in UCINET ("dat") format. Get the node names:
fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
      skip = 4,               # Must look at the data structure to be sure
      #  of this
      header = FALSE,
      nrows = 14)$V1 -> employees

# Figure 2.2(a)
as.matrix(fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
                skip = 41,    # Must look at the data structure to be sure 
                #  of this
                nrows = 14)) -> games_mat
employees -> rownames(games_mat) -> colnames(games_mat)
network(games_mat, directed = FALSE) -> games_net

set.seed(42)
plot.network(games_net,
             label = network.vertex.names(games_net),
             edge.lwd = 2,
             vertex.col = 8,
             vertex.sides = 4,    # This doesn't round the corners of the
             #   verices
             vertex.cex = 3,
             vertex.border = 1,
             pad = 1)             # My plots clip labels otherwise

# Figure 2.2(b)
as.matrix(fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
                skip = 55,    # Must look at the data structure to be sure 
                #  of this
                nrows = 14)) -> friends_mat
employees -> rownames(friends_mat) -> colnames(friends_mat)
network(friends_mat, directed = FALSE) -> friends_net

set.seed(42)
plot.network(friends_net,
             label = network.vertex.names(friends_net),
             edge.lwd = 2,
             vertex.col = 8,
             vertex.sides = 4,    # This doesn't round the corners of the
             #  verices
             vertex.cex = 3,
             vertex.border = 1,
             pad = 1)             # My plots clip labels otherwise


# Figure 2.3
# Read these from the UCI data:
read.ucinet.header("Data/campnet") -> camp_hdr
read.ucinet("Data/campnet") -> camp_mat
camp_hdr$dim.labels[[1]] -> rownames(camp_mat) -> colnames(camp_mat)
network(camp_mat,
        directed = TRUE) -> camp_net

# Get and set the vertex attributes
read.ucinet.header("Data/campsex") -> camp_sex_hdr
read.ucinet("Data/campsex") -> camp_sex_mat
# These got imported as a matrix:
network::set.vertex.attribute(camp_net, "Sex", value = camp_sex_mat[,1]) 
# These got imported as a matrix:
network::set.vertex.attribute(camp_net, "Role", value = camp_sex_mat[,6])

set.seed(42)
plot.network(camp_net,
             label = network.vertex.names(camp_net),
             # Colors rather than shapes to distinguish roles:
             vertex.col = network::get.vertex.attribute(camp_net, "Role") + 1)

# Figure 2.4
# I think I've got the correct data set - the Sampson Monastery data. 
#  However, one needs to be careful, as these data only have 18 of the full
#  set of 25 members. Furthermore, there are 10 different networks of these
#  18, depending upon which question is plotted. So, this may be way off.
# The full data are in Data/sampson.RData.
#
fread("http://moreno.ss.uci.edu/sampson.dat",
      skip = 4,               # Must look at the data structure to be
      #  sure of this
      header = FALSE,
      nrows = 18)$V1 -> monks

# Figure 2.2(a)
as.matrix(fread("http://moreno.ss.uci.edu/sampson.dat",
                skip = 53,    # Must look at the data structure to be sure
                #  of this
                nrows = 18)) -> monks_mat
monks -> rownames(monks_mat) -> colnames(monks_mat)
network(monks_mat, directed = TRUE) -> monks_net
network::set.edge.value(monks_net, "Weight", monks_mat)

set.seed(42)
plot.network(monks_net,
             label = network.vertex.names(monks_net),
             edge.label = network::get.edge.attribute(monks_net, "Weight"),
             pad = 1)

# Matrix 2.1
friends_mat

# Matrix 2.2
# Sadly, I don't think that statnet has a good way to do distances. So, I
#  convert to igraph and work there.
asIgraph(camp_net) -> camp_gr
distances(camp_gr,
          mode = "out") -> camp_dist
V(camp_gr)$vertex.names -> rownames(camp_dist) -> colnames(camp_dist)
camp_dist

# Matrix 2.3
load("Data/DeepSouth.RData")
as.matrix(davisDyn)

## Chapter 3 ####

# Figure 3.2
matrix(data = c(0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 1, 0, 0, 0, 0, 0,
                0, 1, 0, 1, 0, 0, 0, 0,
                0, 1, 1, 0, 1, 0, 0, 0,
                0, 0, 0, 1, 0, 1, 0, 1,
                0, 0, 0, 0, 1, 0, 1, 1,
                0, 0, 0, 0, 0, 1, 0, 1,
                0, 0, 0, 0, 1, 0, 1, 0),
       nrow = 8,
       ncol = 8,
       byrow = TRUE) -> F3.2_mat
network(F3.2_mat,
        directed = FALSE) -> F3.2_net
plot(F3.2_net,
     label = network.vertex.names(F3.2_net),
     vertex.col = 8,
     vertex.sides = 4,
     main = "Tie absent")
sna::betweenness(F3.2_net,
                 gmode = "graph",
                 rescale = TRUE)

add.edge(F3.2_net, tail = 3, head = 8) -> F3.2b_net
plot(F3.2b_net,
     label = network.vertex.names(F3.2b_net),
     vertex.col = 8,
     vertex.sides = 4,
     main = "Tie present")
# I believe that the table of Figure 3.2 is incorrect. Each of the values
#  listed there are only 2/3 of the value they should be. Multiplying the next
#  command by 2/3 yields the values in Figure 3.2. (They are incorrect in the
#  figure because they do not sum to 1.0)
sna::betweenness(F3.2b_net,
                 gmode = "graph",
                 rescale = TRUE)


## Chapter 4 ####

# Figure 4.5
# Yeah, I'm skipping this one because I don't know a quick way to get which
#  subset of IMDB is being shown - the original web-site that maintained those
#  data is defunct - and because I'm too lazy to put in 40 nodes and 96 edges
#  by hand. We'll be back to stuff like this when we look at plotting networks.


## Chapter 5 ####
# This is the chapter were UCINET really comes into play. Hence, there's a lot
#  of stuff here. Essentially, this is the code to do everything that Borgatti
#  et al. do in Chapter 5.

# Matrix format

# Figure 5.1
# It's a MS Excel spreadsheet...you don't need me for this one.

# Figure 5.2
# UCINET specific; R uses a data frame for this

# Figure 5.3 & Matrix 5.1
list("Bill Smith" = c("Carrie Jones", "Doug Johnson", "Eric Morrison"),
     "Eric Morrison" = "Finn Cobb",
     "Doug Johnson" = c("Finn Cobb", "Eric Morrison"),
     "Carrie Jones" = "Finn Cobb") -> nodelist

# I'd bet that there's a good purrr way to do this, but I'm trying
#  to write this script quickly, not elegantly.
sort(c(names(nodelist), "Finn Cobb")) -> nodeNames
matrix(0,
       ncol = length(nodeNames),
       nrow = length(nodeNames),
       dimnames = list(nodeNames,
                       nodeNames)) -> node_mat
for(index in 1:length(nodelist)) {
  for(col in nodelist[[index]]) {
    1 -> node_mat[names(nodelist)[index], col]
    1 -> node_mat[col, names(nodelist)[index]]
  }
}

# Figure 5.4
# Let's do some igraph, just for a change
graph_from_adjacency_matrix(node_mat,
                            mode = "undirected") -> node_gr
set.seed(42)
# In the next, vertex.label.dist, vertex.label.degree, and
#  vertex.size were set by trial and error. (I've done this
#  before, so it didn't take many trials or much error.)
plot(node_gr,
     layout = layout_in_circle(node_gr),
     vertex.label.dist = c(6, 7, 5, 7, 4),
     vertex.shape = "csquare",
     vertex.color = "gray",
     vertex.size = 8,
     vertex.label.degree = c(0, 0, -pi/6, pi, 3*pi/4))

# Figure 5.5
# UCINET specific

# Figure 5.6
# UCINET specific

# Figure 5.7
# Edgelists (without attributes)
matrix(c("Bill Smith", "Eric Morrison",
         "Bill Smith", "Doug Johnson",
         "Bill Smith", "Carrie Jones",
         "Eric Morrison", "Bill Smith",
         "Eric Morrison","Doug Johnson",
         "Eric Morrison", "Finn Cobb",
         "Doug Johnson", "Bill Smith",
         "Doug Johnson", "Eric Morrison",
         "Doug Johnson", "Finn Cobb",
         "Carrie Jones", "Bill Smith",
         "Carrie Jones", "Finn Cobb",
         "Finn Cobb", "Eric Morrison",
         "Finn Cobb", "Doug Johnson",
         "Finn Cobb", "Carrie Jones"),
       ncol = 2,
       nrow = 14,
       byrow = TRUE) -> F5.7_mat
graph_from_edgelist(F5.7_mat,
                    directed = FALSE) -> F5.7_gr

# Figure 5.8
# Edgelists (with attributes)
# I'm going to use a dataset that is similar to the one BEJ
#  display. Where BEJ use "PADGB" for Padgett & Anselm Business,
#  I'm using "Political"; "M" is for Marriage.
#
# This will require some cleaning up, as the edges may run both
#  ways (as in the case of marriage and some political alliances)
#  or not (as in some business/politics, such as loans). I'm not
#  worrying about that here, though, as nothing simply handles 
#  graphs with some directed and some not-directed.

fread(file = here("Data/Florentine Edges.csv"),
      header = TRUE) -> florentine_edges
graph_from_data_frame(florentine_edges,
                      directed = TRUE) -> flo_gr

# Figure 5.9
# Yep, gml is a thing. igraph::read_graph() will read gml
#  (and other) file formats directly. network::read.paj()
#  will read Pajek files, but I don't think that package
#  has anything that will read gml files.

# Figure 5.10
# I'm not interested in typing in the really big matrix, so
#  let's do a smaller one:
matrix(1:16,
       ncol = 4,
       nrow = 4) -> F5.10_mat
F5.10_mat     # (a) Look at it
t(F5.10_mat)  # (b) Look at the transpose

# Figure 5.11
# Just plot the two, BUT USE THE SAME SEED OR THE SAME LAYOUT
#  FOR BOTH PLOTS.
# I won't do that here.

# Imputing (Section 5.4.2)
# First off, you know how I feel about imputing data
# Second, the way that UCINET imputes data is wrong (because
#  there is no good way to impute networked data; imputing
#  data when the samples are not independent is really, really
#  really hard - infinitely harder than doing a good job of
#  imputing non-networked data).
# Third, this is UCINET specific stuff, so I won't do it here
#
# Matrix 5.2
# Matrix 5.3

# Symmetrizing - some analyses need this
# The way to do it is to work with the matrices:
(F5.10_mat + t(F5.10_mat))/2 -> symmetrized_mat
symmetrized_mat
# In addition, sna::symmetrize() will sorta do this type of thing,
#  but you need to give it a "rule". Look at the help file for more
#  information.

# Dichotomizing
# Make the all binary (i.e., drop the weights)
# The easiest thing to do is just ignore the weights. However, if
#  one works with a matrix, then this works:

# Matrix 5.4
matrix(c(0,2,1,0,2,0,0,2,2,2,2,1,1,1,0,
         2,0,3,3,1,4,2,0,2,1,1,2,0,2,0,
         1,3,0,6,1,2,2,1,2,0,2,2,1,1,0,
         0,3,6,0,2,2,1,0,0,0,4,3,1,0,0,
         2,1,1,2,0,1,1,2,1,1,2,1,1,0,0,
         0,4,2,2,1,0,1,2,2,0,2,0,1,0,0,
         0,2,2,1,1,1,0,1,1,0,1,0,2,1,0,
         2,0,1,0,2,2,1,0,2,1,2,0,2,0,0,
         2,2,2,0,1,2,1,2,0,3,3,0,1,1,0,
         2,1,0,0,1,0,0,1,3,0,3,1,0,1,0,
         2,1,2,4,2,2,1,2,3,3,0,0,1,0,0,
         1,2,2,3,1,0,0,0,0,1,0,0,0,0,0,
         1,0,1,1,1,1,2,2,1,0,1,0,0,1,0,
         1,2,1,0,0,0,1,0,1,1,0,0,1,0,0,
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
       ncol = 15,
       nrow = 15,
       byrow = TRUE) -> Mat5.4
# Here's how I got the above: Kindle for Mac - screen capture -
#  print to PDF - OCR - copy and paste into spreadsheet -
#  save as CSV - open in text editor - copy and paste here -
#  add end of line commas. Probably faster than me typing and 
#  certainly more accurate for 225 numbers and commas.

# Matrix 5.5
# Notice that the cutoff for Matrix 5.5 is weight > 1.

matrix(sapply(Mat5.4, function(x) ifelse(x>1, 1, 0)),
       ncol = 15,
       nrow = 15,
       byrow = FALSE) -> Mat5.5
Mat5.5

# Matrix 5.6
# I don't think y'all have done factor analysis. Hence, we should
#  do that before trying to look at Matrix 5.6.
# Stay tuned...

# Figure 5.12
as.matrix(fread(file = here("Data/pv960.csv"),
                header = TRUE)) -> scientists960_mat
colnames(scientists960_mat) -> rownames(scientists960_mat)

fread(file = here("Data/allattrs960.csv"),
      header = TRUE) -> attributes960_df

# Because most people have a lot more papers authored than
#  they have authored with other people (duh!), we ignore
#  the self-loops ("diag = FALSE").
graph_from_adjacency_matrix(scientists960_mat,
                            mode = 'undirected',
                            diag = FALSE) -> sci960_gr
# In igraph, V()$ will access (or add) a vertex property.
attributes960_df$DeptID -> V(sci960_gr)$Department

set.seed(42)
plot(sci960_gr,
     vertex.color = V(sci960_gr)$Department,
     vertex.size = 3,
     vertex.label = NA)

# Matrix 5.7 & Figure 5.13
# OK, I'm skipping this one because there is a problem that
#  I'm not up for at this point. BEJ give data for the
#  Department ID, but not the Department designation. E.g., 
#  "1" rather than "BHS" or whatever Deparment 1 is. Maybe they
#  are in alphabetical order, but the Dept IDs are: 1, 2, 3, 4, 5
#  6, 7, 10, 20, 30. Ugh! If I get bored at some time, I'll
#  figure it out.
#  

# Figure 5.13
# Given Matrix 5.7, we get this graph, ignoring self-loops, and
#  dichotomizing for link strength > 2:
as.matrix(fread(file = here("Data/Matrix 5.7.csv"),
                header = FALSE)) -> Mat5.7
c("BHS", "CCG", "DCL", "ES", "HEW", "IS", "MS", "SRG", "STAT",
  "TAS") -> dept_names
list(dept_names, dept_names) -> dimnames(Mat5.7)
# 1.4 in the next seems to work; BEJ give no guidance
#  on the choice of dichotomizing parameter
matrix(sapply(Mat5.7, function(x) ifelse(x>1.4, 1, 0)),
       ncol = ncol(Mat5.7),
       nrow = nrow(Mat5.7),
       byrow = FALSE) -> Mat5.7_dichot
list(dept_names, dept_names) -> dimnames(Mat5.7_dichot)
graph_from_adjacency_matrix(Mat5.7_dichot,
                            diag = FALSE,
                            mode = "undirected") -> F5.13_gr
set.seed(20200414)
# Won't do anything fancy here.
plot(F5.13_gr,
     vertex.shape = "square",
     vertex.size = 5,
     vertex.label.dist = 1)

# Notice that Section 5.7 doesn't really apply to R; you can do
#  it via the node number as an identifier.
#  Hence, Matrix 5.8, 5.9, and Figure 5.14 we'll skip.


# Figure 5.15
LETTERS[1:6] -> mat_names
c(1,2,2,1,1,2) -> genders
matrix(0,
       ncol = length(genders), 
       nrow = length(genders)) -> gender_mat
list(mat_names, mat_names) -> dimnames(gender_mat)
for(i in 1:length(genders)) {
  for(j in i:length(genders)) {
    if(genders[i] == genders[j]) {
      1 -> gender_mat[i,j] -> gender_mat[j,i]
    }
  }
}
gender_mat

# Figure 5.16
LETTERS[1:6] -> mat_names
c(14, 67, 34, 33, 56, 45) -> ages
matrix(0,
       ncol = length(ages), 
       nrow = length(ages)) -> age_diff_mat
list(mat_names, mat_names) -> dimnames(age_diff_mat)
for(i in 1:length(ages)) {
  for(j in i:length(ages)) {
    abs(ages[i]-ages[j]) -> 
      age_diff_mat[i,j] -> 
      age_diff_mat[j,i]
  }
}
age_diff_mat

# Figure 5.17
LETTERS[1:6] -> mat_names
c(6,10, 3, 5, 9, 4) -> status
matrix(0,
       ncol = length(status), 
       nrow = length(status)) -> status_diff_mat
list(mat_names, mat_names) -> dimnames(status_diff_mat)
for(i in 1:length(status)) {
  for(j in i:length(status)) {
    status[i]-status[j] -> 
      status_diff_mat[i,j] -> 
      status_diff_mat[j,i]
  }
}
status_diff_mat

# Figure 5.18
# Same status and names as Figure 5.17
matrix(0,
       ncol = length(status), 
       nrow = length(status)) -> status_alter_mat
list(mat_names, mat_names) -> dimnames(status_alter_mat)
for(i in 1:length(status)) {
  for(j in 1:length(status)) {
    status[j] -> status_alter_mat[i,j]
  }
}
status_alter_mat

## Still to do: Figure 5.7
## Chapter 6 ####

# Matrix 6.1
matrix(c(0, 206, 429, 1504, 963, 2976, 3095, 2979, 1949,
         206, 0, 233, 1308, 802, 2815, 2934, 2786, 1771,
         419, 233, 0, 1075, 671, 2684, 2799, 2631, 1616,
         1504, 1308, 1075, 0, 1329, 3273, 3053, 2687, 2037,
         963, 802, 671, 1329, 0, 2013, 2142, 2054, 996,
         2976, 2915, 2684, 3273, 2013, 0, 808, 1131, 1307,
         3095, 2934, 2799, 3053, 2142, 808, 0, 379, 1235,
         2979, 2786, 2631, 2687, 2054, 1131, 379, 0, 1059,
         1949, 1771, 1816, 2037, 996, 1307, 1235, 1059, 0
),
nrow = 9,
ncol = 9,
byrow = TRUE) -> city_mat
c("Boston", "NY", "DC", "Miami", "Chicago", "Seattle", "SF", "LA",
  "Denver") -> citynames -> rownames(city_mat) -> colnames(city_mat)

# Figure 6.1
# The "as.dist()" or "dist()" functions change the matrix
#  into something that cmdscale() can work with
as.dist(city_mat) -> city_dist

fit <- cmdscale(city_dist, eig = TRUE, k = 2)
x <- -fit$points[, 1]  # Otherwise, the plot is reversed
y <- fit$points[, 2]   # This should be reversed, but isn't in
#  the text

{
  plot(x, y, pch = 19, xlim = range(x) + c(0, 600))
  text(x, y, pos = 4, labels = citynames)
}

# Table 6.1 & Figure 6.2
# The next works, with warnings
read.ucinet.header("Data/doctorates") -> doc_hdr
read.ucinet("Data/doctorates") -> doc_mat
# A bit of clean-up; remove attributes
NULL -> attr(doc_mat, "title")
NULL -> attr(doc_mat, "date")
doc_hdr$dim.labels[[2]] -> rownames(doc_mat)
doc_hdr$dim.labels[[1]] -> colnames(doc_mat)
doc_mat         # Table 6.1

# This is upside down from BEJ.
ca(doc_mat,
   nd = 2) -> doctorates_ca
plot(doctorates_ca)     # Figure 6.2

# Figure 6.3
matrix(c("a", "b",
         "a", "f",
         "a", "d",
         "b", "d",
         "b", "e",
         "b", "g",
         "c", "d",
         "c", "h",
         "d", "g",
         "d", "f",
         "d", "e",
         "e", "g",
         "f", "g",
         "f", "h",
         "g", "h",
         "h", "r",
         "h", "i",
         "i", "j",
         "i", "n",
         "i", "p",
         "j", "k",
         "j", "l",
         "j", "m",
         "j", "n",
         "j", "p",
         "j", "o",
         "k", "l",
         "k", "p",
         "l", "m",
         "m", "q",
         "n", "o",
         "q", "s"),
       ncol = 2,
       byrow = TRUE) -> F6.3_edges
network(F6.3_edges, directed = FALSE) -> F6.3_net
set.seed(42)
plot.network(F6.3_net,
             label = network.vertex.names(F6.3_net))

# Table 6.2
# Late chapters will talk in more detail about centralities
#  (It is this habit of forward references than made me assign
#  the extra pre-course reading.) Note, also, that BEJ normalize
#  their valued differently than others; hence the multiplications
#  and divisions. (And still, we're off by a bit. I think the
#  issue is some missing edges somewhere...the correspondence
#  analysis will have some differences...oh, well!)
as.data.frame(F6.3_edges) -> F6.3_edgelist
sna::betweenness(F6.3_net) * 
  nrow(F6.3_edgelist) / 100 -> F6.3_between
100 * sna::closeness(F6.3_net) -> F6.3_close
sna::degree(F6.3_net) -> F6.3_degree
# I have no idea about the normalization on the next one, but
#  141 is the correct number
141 * sna::evcent(F6.3_net) -> F6.3_eigen

# Figure 6.3
cbind(F6.3_degree, 
      F6.3_close, 
      F6.3_between, 
      F6.3_eigen) -> F6.3_mat
letters[1:19] -> rownames(F6.3_mat)
c("Degree", "Closeness", "Betweenness", "Eigenvector") -> 
  colnames(F6.3_mat)
F6.3_mat

# Figure 6.4
ca(F6.3_mat, nd = 2) -> F6.3_ca
plot(F6.3_ca)

# Here's one way to put this in a graph:
data.frame("Name" = letters[1:19],
           "Betweenness" = F6.3_between,
           "Closeness" = F6.3_close,
           "Degree" = F6.3_degree,
           "Eigenvector" = F6.3_eigen,
           row.names = NULL) -> F6.3_vert

graph_from_data_frame(d = F6.3_edgelist,
                      directed = FALSE,
                      vertices = F6.3_vert) -> F6.3_gr


# Matrix 6.2 is simply city_mat above.

# Matrices 6.3, 6.4, 6.5 and 6.6 won't be redone here. However, here is the
#  clustering output of Figure 6.5:

hclust(city_dist, method="single") -> hc_single  # The BEJ way
plot(hc_single)
hc_single$height    # Read the dendogram from the bottom to see these in
#  context.

hc_single$labels
hc_single$merge
##      [,1] [,2]
## [1,]   -1   -2   # Merge the first and second labels
## [2,]   -3    1   # Merge the third label into the cluster created in step 1
## [3,]   -7   -8   # Merge the 7th and 8th labels into a cluster
## [4,]   -5    2   # Merge the 5th label into the cluster created in step 2
## [5,]   -6    3   # Merge the 6th label into the cluster created in step 3
## [6,]   -9    4   # Merge the 9th label into the cluster created in step 4
## [7,]    5    6   # Merge the cluster created in step 5 into the cluster
#  created in step 6
## [8,]   -4    7   # Merge the 4th label into the cluster created in step 7

# My preferred method, in general, if I don't know anything else:
hclust(city_dist, method="average") -> hc_ave
# Compare to plot(hc_single); Miami & Denver are better, no?:
plot(hc_ave)
hc_ave$height


## Chapter 7 ####

# Figure 7.1
# There is a random element, so I won't get exactly the same. And...I don't
#  know where the I3 links, for example, come from in BEJ, as I don't see them
#  in the games matrix. Also, I3 doesn't have any links in Figure 7.2 in BEJ.
#  Maybe I3 is just unfortunately placed on a link, but isn't actually
#  connected there. Anyway, I hope you get the idea.
fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
      skip = 4,               # Must look at the data structure to be sure
      #  of this
      header = FALSE,
      nrows = 14)$V1 -> employees
as.matrix(fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
                skip = 41,    # Must look at the data structure to be sure of
                #  this
                nrows = 14)) -> games_mat
employees -> rownames(games_mat) -> colnames(games_mat)
network(games_mat, directed = FALSE) -> games_net

set.seed(42)
plot.network(games_net,
             label = network.vertex.names(games_net),
             mode = "fruchtermanreingold",
             layout.par = list(niter = 1),         # Random start with no 
             #  updating
             edge.lwd = 2,
             vertex.col = 8,
             vertex.sides = 4,    # This doesn't round the corners of the 
             #  vertices
             vertex.cex = 3,
             vertex.border = 1,
             pad = 1)             # My plots clip labels otherwise

#Figure 7.2
# This type of separation can also be accomplished through
#  tcl/tk plotting, but (a) that is interactive (and so it
#  WILL NOT KNIT) and (b) is is a major pain.
set.seed(42)
plot.network(games_net,
             label = network.vertex.names(games_net),
             mode = "fruchtermanreingold",           # Nodes "push"
             #  each other away
             layout.par = list(niter = 500),         # Random start with 500
             #  iterations
             edge.lwd = 2,
             vertex.col = 8,
             vertex.sides = 4,    # This doesn't round the corners of the 
             #  verices
             vertex.cex = 3,
             vertex.border = 1,
             pad = 1)             # My plots clip labels otherwise

# Figure 7.3
fread(file = here("Data/Trade_minerals.csv"),
      header = TRUE) -> trade_df
trade_df[,-1] -> trade_df
as.matrix(trade_df) -> trade_mat
colnames(trade_df) -> rownames(trade_mat)
network(trade_mat, directed = FALSE) -> trade_net

fread(file = here("Data/Trade_Attribute.csv"),
      header = TRUE) -> trade_attr_df
trade_attr_df$SCHOOLS -> x
trade_attr_df$ENERGY -> y
# For some reason, I can't get the coord parameter in plot.network() to work
#  properly. Hence, my work-around is to normalize the x & y variables
x/max(x) -> x_norm
y/max(y) -> y_norm


plot.network(trade_net,
             label = colnames(trade_mat),
             coord = as.matrix(cbind(x_norm, y_norm)),
             # The next line doesn't seem to work as I think it should, making the previous
             #  line & work necessary:
             #             coord = as.matrix(cbind(x, y)),
             label.cex = 0.5)

# Figure 7.4
as.dist(trade_mat) -> trade_dist

fit <- cmdscale(trade_dist, eig = TRUE, k = 2)
x <- -fit$points[, 1]  # Otherwise, the plot is reversed
y <- fit$points[, 2]   # This should be reversed, but isn't in
#  the text

# BEJ obviously must have used different (I suspect additional) data than the
#  data that were available on the book's web-site. Oh, well. I'm not going to
#  investigate that problem right now; this is the correct process.
plot.network(trade_net,
             label = colnames(trade_mat),
             coord = as.matrix(cbind(x, y)),
             label.cex = 0.5)


# Figure 7.5
set.seed(42)
# There is some randomness in all this, so I won't reproduce Figure 7.5
#  perfectly:
plot.network(trade_net,
             label = colnames(trade_mat),
             label.cex = 0.5)

# Figure 7.6
read.ucinet.header("Data/campnet") -> camp_hdr
read.ucinet("Data/campnet") -> camp_mat
camp_hdr$dim.labels[[1]] -> rownames(camp_mat) -> colnames(camp_mat)
network(camp_mat,
        directed = TRUE) -> camp_net

# Get and set the vertex attributes
read.ucinet.header("Data/campsex") -> camp_sex_hdr
read.ucinet("Data/campsex") -> camp_sex_mat
# The next two got imported as a matrix:
network::set.vertex.attribute(camp_net, "Sex", value = camp_sex_mat[,1])  
network::set.vertex.attribute(camp_net, "Betweenness", value = sna::betweenness(camp_net))

# Set edge attributes
(outer(camp_sex_mat[,1], camp_sex_mat[,1])) %% 2 + 1 -> edge_type

c(4, 20) -> shapes       # Squares and approximate circles
set.seed(42)
plot.network(camp_net,
             label = network.vertex.names(camp_net),
             vertex.sides = shapes[network::get.vertex.attribute(camp_net, "Sex")],
             vertex.cex = sqrt(network::get.vertex.attribute(camp_net, "Betweenness")),
             #             edge.lty = edge_type[as.edgelist(camp_net)[,]])
             edge.lty = rep(c(1:2), 12))


# Figure 7.7
read.ucinet.header("Data/campsex") -> camp_sex_hdr
read.ucinet("Data/campsex") -> camp_sex_mat
# The next two got imported as a matrix:
network::set.vertex.attribute(camp_net, "Sex", value = camp_sex_mat[,1])
network::set.vertex.attribute(camp_net, "Betweenness", value = sna::betweenness(camp_net))

# Induced subgraph:
conflict_prefer("%s%", "network")  # Also appears in igraph
camp_net %s% which(network::get.vertex.attribute(camp_net, "Sex") == 1) -> camp_red
set.seed(42)
plot.network(camp_red,
             label = network.vertex.names(camp_red),
             vertex.cex = sqrt(network::get.vertex.attribute(camp_net, "Betweenness")),
             pad = 1)

# I'll skip 7.8 and 7.9...nothing much to learn there

# Figure 7.10
data(karate)
asNetwork(karate) -> karate_net
get.neighborhood(karate_net, 1) -> MrHi_neighbors
c(MrHi_neighbors, 1) -> MrHi_ego
karate_net %s% MrHi_ego -> MrHi_net
plot(MrHi_net,
     label = network.vertex.names(MrHi_net))

# Figure 7.11
# Note that the attribute used here to define shapes disagrees
#  with the plot in BEJ
data(karate)
asNetwork(karate) -> karate_net
get.neighborhood(karate_net, 1) -> MrHi_neighbors
karate_net %s% MrHi_neighbors -> MrHi_net
c(4,20) -> shapes
plot.network(MrHi_net,
             vertex.sides = shapes[network::get.vertex.attribute(MrHi_net, "Faction")],
             vertex.rot = 45,
             label = network.vertex.names(MrHi_net),
             pad = 1)

# Figure 7.12
fread(file = here("Data/davis.csv"),
      header = TRUE) -> davis
as.matrix(davis[,-1]) -> davis_mat
davis$V1 -> rownames(davis_mat)
davis_mat %*% t(davis_mat) -> davis_ord    # This is a way to make an 
#  ordination plot on the rows


graph_from_adjacency_matrix(davis_ord,
                            mode = "undirected",
                            diag = FALSE,
                            weighted = TRUE) -> davis_gr

layout_with_kk(davis_gr,
               weights = 1/E(davis_gr)$weight) -> xy
set.seed(42)
plot(davis_gr,
     layout = xy,
     vertex.size = 6,
     vertex.label.dist = 1)

# Figure 7.13
fread(file = here("Data/davis.csv"),
      header = TRUE) -> davis
as.matrix(davis[,-1]) -> davis_mat
davis$V1 -> rownames(davis_mat)
davis_mat %*% t(davis_mat) -> davis_ord    # This is a way to make an 
#  ordination plot on the rows


graph_from_adjacency_matrix(davis_ord,
                            mode = "undirected",
                            diag = FALSE,
                            weighted = TRUE) -> davis_gr
subgraph.edges(davis_gr, 
               eids = which(E(davis_gr)$weight > 2),
               delete.vertices = FALSE) -> davis_red  # Keeps the isolated 
#  vertices

layout_with_kk(davis_red,
               weights = 1/E(davis_red)$weight) -> xy
set.seed(42)
plot(davis_red,              # Ugly plot
     layout = xy,
     vertex.size = 6,
     vertex.label.dist = 1)

# Figure 7.14
fread(file = here("Data/davis.csv"),
      header = TRUE) -> davis
as.matrix(davis[,-1]) -> davis_mat
davis$V1 -> rownames(davis_mat)
davis_mat %*% t(davis_mat) -> davis_ord    # This is a way to make an 
#  ordination plot on the rows


graph_from_adjacency_matrix(davis_ord,
                            mode = "undirected",
                            diag = FALSE,
                            weighted = TRUE) -> davis_gr

layout_with_kk(davis_gr,
               weights = 1/E(davis_gr)$weight) -> xy
set.seed(42)
plot(davis_gr,
     layout = xy,
     edge.width = E(davis_gr)$weight,
     vertex.size = 6,
     vertex.label.dist = 1)


# Figure 7.15
fread(file = here("Data/Salmon.csv"),
      header = TRUE) -> salmon
as.matrix(salmon[,-1]) -> salmon_mat
salmon$V1 -> rownames(salmon_mat) -> colnames(salmon_mat)

network(salmon_mat, 
        directed = TRUE) -> salmon_net

network::set.edge.attribute(salmon_net, 
                            attrname = "weight", 
                            value = salmon_mat[as.edgelist(salmon_net)[,]])

# In the next, I use the sqrt() to take out some of the variation, and divide
#  by 20 to scale it a little bit better in relation to the size of the 
#  graph.
plot.network(salmon_net,
             arrowhead.cex = sqrt(network::get.edge.attribute(salmon_net, "weight"))/20)

# Figure 7.16
# Figure 7.16(a)
fread(file = here("Data/davis.csv"),
      header = TRUE) -> davis
as.matrix(davis[,-1]) -> davis_mat
davis$V1 -> rownames(davis_mat)
davis_mat %*% t(davis_mat) -> davis_ord    # This is a way to make an 
#  ordination plot on the rows


graph_from_adjacency_matrix(davis_ord,
                            mode = "undirected",
                            diag = FALSE,
                            weighted = TRUE) -> davis_gr
subgraph.edges(graph = davis_gr,
               eids = E(davis_gr)[which(E(davis_gr)$weight > 1)]) -> 
  davis_gr_1


layout_with_kk(davis_gr_1,
               weights = 1/E(davis_gr_1)$weight) -> xy
set.seed(42)
plot(davis_gr_1,
     layout = xy,
     edge.width = E(davis_gr_1)$weight,
     vertex.size = 6,
     vertex.label.dist = 1)


# Figure 7.16(b)
fread(file = here("Data/davis.csv"),
      header = TRUE) -> davis
as.matrix(davis[,-1]) -> davis_mat
davis$V1 -> rownames(davis_mat)
davis_mat %*% t(davis_mat) -> davis_ord    # This is a way to make an 
#  ordination plot on the rows


graph_from_adjacency_matrix(davis_ord,
                            mode = "undirected",
                            diag = FALSE,
                            weighted = TRUE) -> davis_gr
subgraph.edges(graph = davis_gr,
               eids = E(davis_gr)[which(E(davis_gr)$weight > 2)]) -> 
  davis_gr_2


layout_with_kk(davis_gr_2,
               weights = 1/E(davis_gr_2)$weight) -> xy
set.seed(42)
plot(davis_gr_2,
     layout = xy,
     edge.width = E(davis_gr_2)$weight,
     vertex.size = 6,
     vertex.label.dist = 1)


# Figure 7.16(c)
fread(file = here("Data/davis.csv"),
      header = TRUE) -> davis
as.matrix(davis[,-1]) -> davis_mat
davis$V1 -> rownames(davis_mat)
davis_mat %*% t(davis_mat) -> davis_ord    # This is a way to make an 
#  ordination plot on the rows


graph_from_adjacency_matrix(davis_ord,
                            mode = "undirected",
                            diag = FALSE,
                            weighted = TRUE) -> davis_gr
subgraph.edges(graph = davis_gr,
               eids = E(davis_gr)[which(E(davis_gr)$weight > 3)]) -> 
  davis_gr_3


layout_with_kk(davis_gr_3,
               weights = 1/E(davis_gr_3)$weight) -> xy
set.seed(42)
plot(davis_gr_3,
     layout = xy,
     edge.width = E(davis_gr_3)$weight,
     vertex.size = 6,
     vertex.label.dist = 1)


# Figure 7.16(d)
fread(file = here("Data/davis.csv"),
      header = TRUE) -> davis
as.matrix(davis[,-1]) -> davis_mat
davis$V1 -> rownames(davis_mat)
davis_mat %*% t(davis_mat) -> davis_ord    # This is a way to make an 
#  ordination plot on the rows


graph_from_adjacency_matrix(davis_ord,
                            mode = "undirected",
                            diag = FALSE,
                            weighted = TRUE) -> davis_gr
subgraph.edges(graph = davis_gr,
               eids = E(davis_gr)[which(E(davis_gr)$weight > 4)]) -> 
  davis_gr_4


layout_with_kk(davis_gr_4,
               weights = 1/E(davis_gr_4)$weight) -> xy
set.seed(42)
plot(davis_gr_4,
     layout = xy,
     edge.width = E(davis_gr_4)$weight,
     vertex.size = 6,
     vertex.label.dist = 1)

# Figure 7.17
fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
      skip = 4,               # Must look at the data structure to be sure
      #  of this
      header = FALSE,
      nrows = 14)$V1 -> employees
as.matrix(fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
                skip = 41,    # Must look at the data structure to be sure of
                #  this
                nrows = 14)) -> games_mat
employees -> rownames(games_mat) -> colnames(games_mat)
network(games_mat, directed = FALSE) -> games_net

set.seed(42)
plot.network(games_net,
             label = network.vertex.names(games_net),
             mode = "fruchtermanreingold",           # Nodes "push"
             #  each other away
             layout.par = list(niter = 500),         # Random start with 500
             #  iterations
             edge.lwd = 2,
             vertex.col = 8,
             vertex.sides = 4,    # This doesn't round the corners of the 
             #  verices
             vertex.cex = 3,
             vertex.border = 1,
             pad = 1)             # My plots clip labels otherwise

# Figure 7.18
fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
      skip = 4,               # Must look at the data structure to be sure
      #  of this
      header = FALSE,
      nrows = 14)$V1 -> employees
as.matrix(fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat",
                skip = 55,    # Must look at the data structure to be sure of
                #  this
                nrows = 14)) -> conflict_mat
employees -> rownames(conflict_mat) -> colnames(conflict_mat)
network(conflict_mat, directed = FALSE) -> conflict_net

set.seed(42)
plot.network(conflict_net,
             label = network.vertex.names(conflict_net),
             mode = "fruchtermanreingold",           # Nodes "push"
             #  each other away
             layout.par = list(niter = 500),         # Random start with 500
             #  iterations
             edge.lwd = 2,
             vertex.col = 8,
             vertex.sides = 4,    # This doesn't round the corners of the 
             #  verices
             vertex.cex = 3,
             vertex.border = 1,
             pad = 1)             # My plots clip labels otherwise

# Matrix 7.1
matrix(c(1.000, 0.684, 0.483, 0.440, 0.300,
         0.684, 1.000, 0.582, 0.543, 0.335,
         0.483, 0.582, 1.000, 0.613, 0.341,
         0.440, 0.543, 0.613, 1.000, 0.371,
         0.300, 0.335, 0.341, 0.371, 1.000),
       nrow = 5,
       ncol = 5,
       byrow = TRUE,
       dimnames = list(c("T1", "T2", "T3", "T4", "T5"),
                       c("T1", "T2", "T3", "T4", "T5"))) -> BandB_mat
BandB_diss = 1 - BandB_mat    # Dissimilarity (not similarity!) matrix

# Figure 7.19
cmdscale(as.dist(BandB_diss), eig = TRUE, k = 2) -> BandB_fit

x <- BandB_fit$points[, 1]
y <- BandB_fit$points[, 2]
# Don't really need a network for this, eh?
{
  plot(x, y, pch = 0,
       ylim = c(-0.3, 0.3),     # Set the ylim a bit big so the label fits
       ann = FALSE,             # Do not annotate the axes
       axes = FALSE,            # Do not show the axes
       frame.plot = TRUE)
  text(x, y,
       labels = colnames(BandB_mat),
       pos = 3)
}

# Figure 7.20 (almost)
# I'm not convinced that these data are intepreted correctly.
#  Even when I dichotomize them, they still don't come anywhere near
#  what BEJ have. Oh, well.
read.xlsx(here("Data/burkhardt.xlsx"),
          sheetName = "T1") -> burkhardt_t1
burkhardt_t1[,1] -> rownames(burkhardt_t1)
burkhardt_t1[,-1] -> burkhardt_t1
as.matrix(burkhardt_t1) -> burkhardt_t1

network(burkhardt_t1,
        directed = FALSE) -> burk_t1_net
set.seed(42)
plot(burk_t1_net) # Skip the label of R53

# Figure 7.21
# I'm not convinced that these data are intepreted correctly.
#  Even when I dichotomize them, they still don't come anywhere near
#  what BEJ have. Oh, well.
read.xlsx(here("Data/burkhardt.xlsx"),
          sheetName = "T5") -> burkhardt_t5
burkhardt_t5[,1] -> rownames(burkhardt_t5)
burkhardt_t5[,-1] -> burkhardt_t5
as.matrix(burkhardt_t5) -> burkhardt_t5

network(burkhardt_t5,
        directed = FALSE) -> burk_t5_net
set.seed(42)
plot(burk_t5_net) # Skip the label of R53

# Figure 7.23
# There is a problem in the dataset (at least as far
#  as ca() is concerned: One row (JEAN2) is all zeroes. This is,
#  mathematically, bad, and ca() doesn't appear to have alternate
#  methods, so there will be matrix problems.
# The only way around this that I can see is to delete JEAN
#  from the analyses. Based on the BEJ plot, they may have 
#  JEAN2 at (0,0), but that shouldn't be used for missing
#  data.
# Hence, I'm going to delete JEAN and JEAN2 from consideration.
fread(file = here("Data/workshop.csv"),
      header = TRUE) -> workshop
workshop[-27,] -> workshop  # JEAN2 row
workshop[-10,] -> workshop  # JEAN row
workshop$V1 -> work_names
as.matrix(workshop[,-1]) -> workshop
workshop[,-10] -> workshop  # JEAN column
work_names -> rownames(workshop)

ca(workshop, nd = 2) -> work_ca
# I find this method easier for this than trying to
#  fiddle with plot.ca().
{
  plot(work_ca$rowcoord[,1],work_ca$rowcoord[,2],
       cex = 0.7,
       ann = FALSE,
       axes = FALSE,
       frame.plot = TRUE)
  text(work_ca$rowcoord[,1],work_ca$rowcoord[,2],
       labels = work_names,
       pos = 4,
       cex = 0.6)
  arrows(x0 = work_ca$rowcoord[1:16,1],
         y0 = work_ca$rowcoord[1:16,2],
         x1 = work_ca$rowcoord[17:32,1],
         y1 = work_ca$rowcoord[17:32,2],
         length = 0.08)  # This is annoying; it depends a lot
  #  on the size of the picture.
  # Zoom in on the final picture.
}

# Figure 7.24
fread(file = here("Data/Supremeall.csv"),
      header = TRUE) -> supreme
supreme$V1 -> sup_names
as.matrix(supreme[,-1]) -> supreme
sup_names -> rownames(supreme)

ca(supreme, nd=2) -> supreme_ca
plot(supreme_ca) # Again, backwards from BEJ

# Figure 7.25
fread(file = here("Data/Supremeall.csv"),
      header = TRUE) -> supreme
supreme$V1 -> sup_names
as.matrix(supreme[,-1]) -> supreme
sup_names -> rownames(supreme)

ca(supreme, nd=2) -> supreme_ca

# ID Renquist
seq(from = 1, to = 90, by = 9) -> renquist
{
  plot(supreme_ca$rowcoord[,1],supreme_ca$rowcoord[,2],
       cex = 0.7,
       ann = FALSE,
       axes = FALSE,
       frame.plot = TRUE)
  text(supreme_ca$rowcoord[,1],supreme_ca$rowcoord[,2],
       labels = sup_names,
       pos = 4,
       cex = 0.6)
  arrows(x0 = supreme_ca$rowcoord[renquist[1:9],1],
         y0 = supreme_ca$rowcoord[renquist[1:9],2],
         x1 = supreme_ca$rowcoord[renquist[2:10],1],
         y1 = supreme_ca$rowcoord[renquist[2:10],2],
         length = 0.08)  # This is annoying; it depends a lot
  #  on the size of the picture.
  # Zoom in on the final picture.
}

# Figure 7.26
# Community structure does this too.
# I won't do this one right now because I can't find the
#  data needed to do this. However...
# The function you want is alphahull::ahull(). Once
#  you have the coordinates of the liberal, swing, and
#  conservative nodes, you can feed that into ahull(),
#  and get out the resulting boundaries.

## Chapter 8 ####
# Figure 8.1
# I already know the file structure...
as.matrix(fread(here("Data/PADGETT_Business.csv"),
                select = 2:17,
                header = TRUE)) -> padgb
fread(here("Data/PADGETT_Business.csv"),
      select = 1) -> row_names
t(row_names) -> rownames(padgb) -> colnames(padgb)

as.matrix(fread(here("Data/PADGETT_Marriage.csv"),
                select = 2:17,
                header = TRUE)) -> padgm
fread(here("Data/PADGETT_Marriage.csv"),
      select = 1) -> row_names
t(row_names) -> rownames(padgm) -> colnames(padgm)

array(dim = c(2, 16, 16)) -> padg_arr
padgm -> padg_arr[1,,]
padgb -> padg_arr[2,,]

set.seed(24322)             # This is different than a UCINET seed of 24322

qaptest(padg_arr,
        FUN = gcor,
        g1 = 1,
        g2 = 2,
        reps = 50000) -> padg_qap   # BEJ do 50k reps
# The next prints out differently than does the BEJ. But still, there
#  is similar information: The non-trivial p-value, and the non-trivial
#  correlation is present.
# qaptest() does only the Pearson correlation; you can worry about the others
#  later.
summary(padg_qap)

# Figures 8.2 & 8.3
# These are a lot different in R. Additionally, it doesn't
#  appear that BEJ have defined "reciprocity" or "transitivity" 
#  yet. (They have sort of defined the former in Chapter 5, but
#  we didn't officially do that chapter!)
# Hence, we'll do this elsewhere.

# Figure 8.4
# Note that the campnet data are already dichotomized to
#  the top three
read.ucinet.header("Data/campnet") -> camp_hdr
read.ucinet("Data/campnet") -> camp_mat
camp_hdr$dim.labels[[1]] -> rownames(camp_mat) -> colnames(camp_mat)
network(camp_mat,
        directed = TRUE) -> camp_net
set.seed(42)
plot.network(camp_net,
             label = network.vertex.names(camp_net),
             pad = 1)

# Matrix 8.1
# Age differences in the Krackhardt data
fread(here("Data/High-Tec-Attributes.csv")) -> hitech_attr
matrix(0,
       nrow = nrow(hitech_attr),
       ncol = nrow(hitech_attr)) -> age_diff_mat
for(i in 1:nrow(hitech_attr)) {
  for(j in 1:nrow(hitech_attr)) {
    hitech_attr$AGE[i] - hitech_attr$AGE[j] ->
      age_diff_mat[i,j]
  }
}

# Figure 8.5
as.matrix(
  fread(here("Data/Krack-High-Tec_advice_Friendship_Reports_To.csv"),
        header = TRUE)) -> reports_mat
reports_mat[,-1] -> reports_mat
array(0,
      dim = list(2, nrow(hitech_attr), nrow(hitech_attr))) -> 
  hi_tech_arr
reports_mat -> hi_tech_arr[1,,]
age_diff_mat -> hi_tech_arr[2,,]
set.seed(42)
qaptest(hi_tech_arr,
        gcor,
        g1 = 1,
        g2 = 2,
        reps = 5000) -> F8.5
summary(F8.5)

# Matrix 8.2
matrix(rep(hitech_attr$AGE,
           nrow(hitech_attr)),
       nrow = nrow(hitech_attr),
       ncol = nrow(hitech_attr),
       byrow = TRUE) -> M8.2
M8.2

# Figure 8.6
read.ucinet.header("Data/campnet") -> camp_hdr
read.ucinet("Data/campnet") -> camp_mat

read.ucinet.header("Data/campsex") -> camp_attr_hdr
read.ucinet("Data/campsex") -> camp_attr_mat
nrow(camp_attr_mat) -> rows

matrix(0, nrow = rows, ncol = rows) -> sex_diff_mat
for(i in 1:rows) {
  for(j in i:rows) {
    if(camp_attr_mat[i,1] == camp_attr_mat[j,1]) {
      1 -> sex_diff_mat[i,j] -> sex_diff_mat[j,i]
    }
  }
}

array(0,
      dim = list(2, rows, rows)) -> camp_arr
camp_mat -> camp_arr[1,,]
sex_diff_mat -> camp_arr[2,,]
qaptest(camp_arr,
        gcor,
        g1 = 1,
        g2 = 2,
        reps = 50000) -> F8.6
summary(F8.6)

# Figure 8.7
# Well, if I knew which data were being used, I'd do this one.
# However, BEJ are really using this example simply to 
#  demonstrate the differences between resampling approaches
#  to p-values and the classical approaches to p-values. I
#  believe that you already know this one, so no loss.

# Sections 8.7 and 8.8
# ERGM and SAOM are handled much differently in R than in UCINET.
#  Hence, trying to just recreate the figures here is not only
#  frustrating for the instructor, it is likely to be of limited
#  help. See the appropriate folders in the Blackboard course site
#  for this material.

## Chapter 9 ####

# Figure 9.1
# Nothing to see here...

# Figure 9.2
# No data

# Table 9.1 & Figure 9.3
# Read these from the UCI data:
read.ucinet.header("Data/campnet") -> camp_hdr
read.ucinet("Data/campnet") -> camp_mat
camp_hdr$dim.labels[[1]] -> rownames(camp_mat) -> colnames(camp_mat)

# Get and set the vertex attributes
read.ucinet.header("Data/campsex") -> camp_sex_hdr
read.ucinet("Data/campsex") -> camp_sex_mat  # Sex is column 1

# We'll construct the table by removing from the matrix
#  the links that don't belong, and then computing the density. 
#  This is best done via the matrices. To make life easier, we should
#  break the matrices into M/F blocks first. By good fortune, the
#  matrices are already sorted into M (camp_sex_mat[,1]==2) & F 
#  (camp_sex_mat[,1] == 1). The first 8 are F, the last 10 are M

camp_mat[1:8, 1:8] -> FF_mat
camp_mat[1:8, 9:18] -> FM_mat
camp_mat[9:18, 9:18] -> MM_mat
camp_mat[9:18, 1:8] -> MF_mat


network(FF_mat, directed = TRUE) -> FF_net
network(FM_mat, directed = TRUE) -> FM_net
network(MM_mat, directed = TRUE) -> MM_net
network(MF_mat, directed = TRUE) -> MF_net

list(FF_net, FM_net, MM_net, MF_net) -> sub_mat
gden(sub_mat)
# Notice: the FF & MM densities are the same as BEJ, but the
#  FM & MF are approximately half of what BEJ report; I haven't
#  tracked down that factor, but it is undoubtedly due to the averaging
#  issues in the definition.

network(camp_mat,
        directed = TRUE) -> camp_net
network::set.vertex.attribute(camp_net, "Sex", camp_sex_mat[,1])
set.seed(42)
plot.network(camp_net,
             label = network.vertex.names(camp_net),
             vertex.sides = 
               ifelse(camp_net %v% "Sex" == 2, 4, 50),
             vertex.rot = 45,   # squares, not diamonds
             vertex.cex = 3,
             pad = 1)

# Figure 9.4
# Nothing to see here

# Figure 9.5
# Nothing to see here. The question is: Do you need to know this
#  nomenclature? Unless you are going to be engaging with the 
#  literature, I suspect not. 

# Table 9.2 & Figure 9.6
# The network data used to create these are not provided. However,
#  sna::triad.census() applied to the appropriate network would give
#  those results. 

# Figure 9.7
# Nothing to see here

# Figure 9.8
# Nothing to see here

# Figures 9.9 & 9.10
# No data provided; nothing to see here

## Chapter 10 ####
# FWIW, package:CINNA does a bunch of centralities.
#  For example, try this (if CINNA is installed):
CINNA::proper_centralities(zachary)
# This indicates that there are 49 appropriate centrality measures
#  for the undirected Zachary karate measruement. For directed
#  networks:
data("foodwebs")    # A list of food web igraphs
CINNA::proper_centralities(foodwebs$ChesLower)
# This indicates that there are 51 appropriate centrality measures
#  for the directed Lower Chesapeake food web
# Obviously, centrality is not a clearly defined thing yet.
#  package:keyplayer does even more stuff.

# Figure 10.1
read.ucinet.header(here("Data/campnet")) -> camp_hdr
read.ucinet(here("Data/campnet")) -> camp_mat

network(camp_mat,
        directed = FALSE) -> camp_net  # directed = FALSE symmetrizes it
(sna::degree(camp_net)) / 2 -> deg     # Don't double count!

# We know the first 8 are F, and the second 8 are M, from previous
#  work.
data.frame("Degre" = deg,
           "nDegr" = deg / sum(deg),
           "Sex" = c(rep("F",8), rep("M",10))) -> degree_df
rownames(camp_mat) -> rownames(degree_df)

# Once again, I don't get the same normalization that BEJ do. Alas
degree_df[,1:2]

# Figure 10.2
# First approach, lmPerm::aovp
summary(aovp(Degre ~ Sex, data = degree_df))  # No difference

# Second approach: Permutation via base R
mean(degree_df$Degre[degree_df$Sex=="M"]) - 
  mean(degree_df$Degre[degree_df$Sex=="M"]) -> obsdiff

resamp_means = function(raw_data, labels){
  resample_labels = sample(labels)
  unique(labels) -> uniq_lab
  resample_diff = mean(raw_data[resample_labels==uniq_lab[1]]) - 
    mean(raw_data[resample_labels==uniq_lab[2]])
  return(resample_diff)
}

1e4 -> N
set.seed(42)
replicate(N,resamp_means(degree_df$Degre,
                         degree_df$Sex)) -> nulldist
{
  hist(nulldist, col="cyan", breaks = 30)
  abline(v = obsdiff, col = "red")   # Looks non-significant
}
# (hist() has seme strange behaviors in terms of those gaps.)

length(which(nulldist > obsdiff)) / N

# Figure 10.3
# Nothing to see here...no data...just an example

# Figure 10.4
# Nothing to see here...no data...just an example

# Figure 10.5
# This really isn't used, but still...
matrix(c("Z6", "Y2",
         "Z5", "Y2",
         "Z4", "Y2",
         "Y2", "X",
         "Y1", "X",
         "Z1", "Y1",
         "Z2", "Y1",
         "Z3", "Y1"),
       ncol = 2,
       nrow = 8,
       byrow = TRUE) -> F10.5_mat
network(F10.5_mat,
        directed = FALSE) -> F10.5_net
set.seed(42)
plot(F10.5_net,
     displaylabels = TRUE)

# Figure 10.6
matrix(c(1,4,
         4,5,
         3,5,
         2,5),
       ncol = 2,
       nrow = 4,
       byrow = TRUE) -> F10.6_mat
network(F10.6_mat,
        directed = TRUE) -> F10.6_net
set.seed(42)
plot(F10.6_net,
     vertex.cex = 2.5,
     vertex.col = "gray",
     displaylabels = TRUE,
     arrowhead.cex = 2)

# Table 10.1
# Some versions of eigenvector centrality allow one to put in
#  a value for _beta_ (really, 1/beta = lambda). However, BEJ's
#  beta centrality is really similar to Bonacich Power Centrality,
#  so we'll use that here.
matrix(c(1,4,
         4,5,
         3,5,
         2,5),
       ncol = 2,
       nrow = 4,
       byrow = TRUE) -> F10.7_mat
graph_from_edgelist(F10.7_mat) -> F10.7_gr
# Because the scaling of power_centrality doesn't match
#  the scaling of BEJ's beta centrality, I'll just do a few to
#  show the similarities. Notice that the ratios of the scores
#  between the nodes is pretty much the same for power_centrality()
#  and for BEJ; just the absolute values are different.
# First, for the "Out scores"
1 -> beta
power_centrality(F10.7_gr)

1000 -> beta
power_centrality(F10.7_gr,
                 exponent = beta)

# Now, let's do the "In scores"
matrix(c(F10.7_mat[,2], F10.7_mat[,1]),
       ncol = 2,
       nrow = 4,
       byrow = FALSE) -> F10.7_in_mat
graph_from_edgelist(F10.7_in_mat) -> F10.7_in_gr

1 -> beta
power_centrality(F10.7_in_gr)

1000 -> beta
power_centrality(F10.7_in_gr,
                 exponent = beta)

# Figure 10.7
# Data not given

## Chapter 11 - still to do ####

# Figure 11.1
matrix(c(1,2,
         1,3,
         1,4,
         2,3,
         2,4,
         3,4,
         3,5,
         4,5,
         5,6,
         6,7,
         7,8,
         7,9,
         7,10,
         8,10,
         9,10),
       ncol = 2,
       nrow = 15,
       byrow = TRUE) -> F11.1_edges
network(F11.1_edges, directed = FALSE) -> F11.1_gr
set.seed(42)
plot(F11.1_gr, displaylabels = TRUE)

# Matrix 11.1
fread(here("Data/WIRING_RDGAM.csv"),
      header = TRUE) -> games_mat
as.matrix(games_mat[,-1]) -> games_mat
colnames(games_mat) -> rownames(games_mat)
games_mat

# Matrix 11.2
# FWIW, here are the 5 cliques that get listed by BEJ (by node index,
#  not node name):
graph_from_adjacency_matrix(games_mat,
                            mode = "directed") -> games_gr

# Can do this with sna::clique.census() too.
max_cliques(games_gr,       # Ignore the warning
            min = 3) -> mc  # BEJ ignore all less than size 3
matrix(0,
       nrow = nrow(games_mat),
       ncol = ncol(games_mat),
       dimnames = list(rownames(games_mat),
                       colnames(games_mat))) -> overlaps

for(index in 1:(length(mc))) {
  for(r in 1:(length(mc[[index]]))) {
    for(c in 1:(length(mc[[index]]))) {
      overlaps[mc[[index]][r],mc[[index]][c]] + 1 -> 
        overlaps[mc[[index]][r],mc[[index]][c]]
    }
  }
}

overlaps

# Figure 11.2
# Start with overlaps

# overlaps is the inverse (or opposite) of the distance, so let's
#  do this:
as.dist(max(overlaps) - overlaps) -> dist_mat
hclust(dist_mat) -> c_map
# Notice that this puts S4 and S2 in, which BEJ's approach 
#  never does. Otherwise, pretty much the same thing, although 
#  the height (levels) of joining are a bit different.
plot(c_map)

# Matrix 11.3
fread(here("Data/WIRING_RDGAM.csv"),
      header = TRUE) -> games_mat
as.matrix(games_mat[,-1]) -> games_mat
colnames(games_mat) -> rownames(games_mat)
colnames(games_mat) -> nam

graph_from_adjacency_matrix(games_mat,
                            mode = "directed") -> games_gr

# Can do this with sna::clique.census() too.
max_cliques(games_gr,       # Ignore the warning
            min = 3) -> mc  # BEJ ignore all less than size 3

# Make the matrix
matrix(0,
       ncol = length(mc),
       nrow = gorder(games_gr),
       dimnames = list(rownames(games_mat), 1:5)) -> bimodal_mat

# Fill the matrix
for(clique in 1:(length(mc))) {
  for(vert in 1:(length(mc[[clique]]))) {
    1 -> bimodal_mat[mc[[clique]][vert], clique]
  }
}

# Now, to fill in the rest, column by column
for(clique in 1:length(mc)) {        # For each clique
  # get the neighbors of the first vertex in the clique:
  neighbors(games_gr, mc[[clique]][1]) -> neigh_vert
  # get the neighbors of each other vertex in the clique, and
  #  add them to the vector:
  for(vert in 2:length(mc[[clique]])) {      
    c(neigh_vert,
      neighbors(games_gr, mc[[clique]][vert])) -> neigh_vert
  }
  # Get the vertices that aren't in the clique
  sort(difference(neigh_vert, mc[[clique]])) -> adds
  # BEJ: "and values in between give the number of ties
  #  which connect the actor to the clique divided by the
  #  number of ties required for the actor to be a clique
  #  member." Notice, however, that BEJ don't actually
  #  follow their rule in the creation of Matrix 11.3. I will
  #  create Matrix 11.3, but the denominator in the next
  #  line is 1 less than is should be.
  unname(table(adds)) / length(mc[[clique]]) -> vals
  # Put the appropriate values in the appropirate positions
  vals -> bimodal_mat[sort(unique(adds)), clique]
}

# Figure 11.3
# bimodal_mat is an affiliation matrix. We can create the network
#  from that. See the code for Matrix 11.3 for the creation of
#  bimodal_mat.
graph_from_incidence_matrix(bimodal_mat) -> bimodal_gr
set.seed(42)
plot(bimodal_gr,
     vertex.shape = ifelse(V(bimodal_gr)$type == TRUE,
                           "square", "circle"),
     vertex.rot = 45,
     vertex.col = c(rep("gray", 19), rep("cyan", 5)),
     vertex.label.dist = 1.5,
     vertex.label.color = "black")

# Figure 11.4
data(karate)
max_cliques(karate,
            min = 3) -> karate_mc

matrix(0,
       nrow = gorder(karate),
       ncol = gorder(karate),
       dimnames = list(V(karate)$name,
                       V(karate)$name)) -> overlaps

for(index in 1:(length(karate_mc))) {
  for(r in 1:(length(karate_mc[[index]]))) {
    for(c in 1:(length(karate_mc[[index]]))) {
      overlaps[karate_mc[[index]][r],
               karate_mc[[index]][c]] + 1 -> 
        overlaps[karate_mc[[index]][r],
                 karate_mc[[index]][c]]
    }
  }
}

# This produces a similar graph to BEJ; that has to do
#  with the next command. (I suspect that if I did something
#  like using rowSums or something else to adjust things,
#  the results would be more like BEJ.)
as.dist(max(overlaps) - overlaps) -> dist_mat
hclust(dist_mat) -> c_map
plot(c_map)

# Figure 11.5

# Figure 11.6

# Figure 11.7
# Nothing to see here...just a picture for the paragraph

# Matrix 11.4
# Everything but the lines in the printout:
matrix(0,
       nrow = 12,
       ncol = 12) -> M11.4
1 -> M11.4[1:3, 1:3]
1 -> M11.4[4:8, 4:8]
1 -> M11.4[9:12, 9:12]
0 -> diag(M11.4)
M11.4

# Figure 11.8
# Girvan-Newman is edge_betweenness in igraph:
data(karate)
cluster_edge_betweenness(karate) -> eb_community # Ignore warning

# This produced 6 communities! BEJ only want only 2. Looking at
#  the membership list, I think that putting 1, 2, 3, and 4 into
#  one community and 5 and 6 in the other will work out just fine.
membership(eb_community) -> comm_memb
1 -> comm_memb[which(comm_memb == 2)]
1 -> comm_memb[which(comm_memb == 3)]
1 -> comm_memb[which(comm_memb == 4)]
2 -> comm_memb[which(comm_memb == 5)]
2 -> comm_memb[which(comm_memb == 6)]


# Figure 11.9 addendum
# Look what is in the igraph documentation. Notice that it has
#  the same degree structure as the UCINET data, but with
#  different names attached. I don't know which is actually
#  correct, but let's not be distracted by that right now; check
#  out the method.

## The science camp network
camp <- graph_from_literal(Harry:Steve:Don:Bert - Harry:Steve:Don:Bert,
                           Pam:Brazey:Carol:Pat - Pam:Brazey:Carol:Pat,
                           Holly   - Carol:Pat:Pam:Jennie:Bill,
                           Bill    - Pauline:Michael:Lee:Holly,
                           Pauline - Bill:Jennie:Ann,
                           Jennie  - Holly:Michael:Lee:Ann:Pauline,
                           Michael - Bill:Jennie:Ann:Lee:John,
                           Ann     - Michael:Jennie:Pauline,
                           Lee     - Michael:Bill:Jennie,
                           Gery    - Pat:Steve:Russ:John,
                           Russ    - Steve:Bert:Gery:John,
                           John    - Gery:Russ:Michael)
campBlocks <- cohesive_blocks(camp)
campBlocks

plot(campBlocks, camp, vertex.label=V(camp)$name, margin=-0.2,
     vertex.shape="rectangle", vertex.size=24, vertex.size2=8,
     mark.border=1, colbar=c(NA, NA,"cyan","orange") )

# Matrix 11.5

# Figure 11.10

# Matrix 11.6

# Figure 11.11

# Matrix 11.7

# Figure 11.12


## Chapter 12 - still to do ####

# Figure 12.1
# Let's just do Relation 2
c(1, 1, 1 + sqrt(3)/2, 1 + sqrt(3)/2, 1 - sqrt(3)/2) -> x
c(0, 1, 1, 0, 0.5) -> y
{
  plot(x, y)
  text(x + c(0.05, 0.05, 0, 0, 0.1),
       y + c(0.1, -0.1, -0.1, 0.1, 0),
       labels = c(1, 2, 3, 4, 5))
  arrows(x0 = c(x[1], x[1], x[1], x[2], x[2], x[2], x[5], x[5]),
         y0 = c(y[1], y[1], y[1], y[2], y[2], y[2], y[5], y[5]),
         x1 = c(x[2], x[3], x[4], x[1], x[3], x[4], x[1], x[2]),
         y1 = c(y[2], y[3], y[4], y[1], y[3], y[4], y[1], y[2]),
         length = 0.1,
         angle = 15,
         ann = FALSE)
}

# Matrix 12.1
matrix(c(0, 1, 0, 1, 1,
         1, 0, 1, 1, 0,
         1, 0, 0, 1, 1,
         1, 0, 1, 0, 0,
         0, 1, 0, 0, 1),
       ncol = 5,
       nrow = 5,
       byrow = TRUE) -> F12.2_mat
as.character(1:5) -> rownames(F12.2_mat) -> colnames(F12.2_mat)
F12.2_mat

# Matrix 12.2
as.matrix(fread("http://moreno.ss.uci.edu/sampson.dat",
                skip = 126,
                nrows = 18,
                header = FALSE)) -> esteem_mat
as.matrix(fread("http://moreno.ss.uci.edu/sampson.dat",
                skip = 144,
                nrows = 18,
                header = FALSE)) -> disesteem_mat
c("ROMULALD", "BONAVENTURE", "AMBROSE", "BERTHOLD", "PETER", 
  "LOUIS", "VICTOR", "WINFRID", "JOHN", "GREG", 
  "HUGH",  "BONIFACE", "MARK", "ALBERT", "AMAND",
  "BASIL", "ELIAS", "SIMPLICIUS") -> row_names
c("RO", "BO", "AM", "BE", "PE", 
  "LO", "VI", "WI", "JO", "GR", "HU", 
  "BO", "MA", "AL", "AM", "BA", "EL", 
  "SI") -> col_names
row_names -> rownames(esteem_mat) -> rownames(disesteem_mat)
col_names -> colnames(esteem_mat) -> colnames(disesteem_mat)

# (a)
esteem_mat

# (b)
disesteem_mat

# Matrix 12.3
# Append the matrices from Matrix 12.2
array(0, dim = list(2, 18, 18)) -> sampson_arr
esteem_mat -> sampson_arr[1,,]
disesteem_mat -> sampson_arr[2,,]

# Sadly, the next computes Matrix 12.3 on the fly, but never
#  makes it available externally. However, we can get the
#  same type of analysis as Matrix 12.3 by looking at the 
#  dendogram in Figure 12.2.
equiv.clust(sampson_arr) -> equiv_res

# Howver, we can calculate the desired distances:
sedist(sampson_arr,
       method = "euclidean",
       joint.analysis = TRUE) -> distances
row_names -> rownames(distances)
col_names -> colnames(distances)
# Notice that the next gives nubmers that are much like
#  those in BEJ in their variation. (They are mostly smaller,
#  but I'm only concerned about ranking, not absolute value.)
round(distances, 1)

# Figure 12.2
plot(equiv_res)
# Looking at the diagram at a heigh of about 25:
# Cluster 1: 8, 14, 15, 18 (Winfrid, Albert, Amand, Simplicius)
# Cluster 2: 11, 12, 13 (High, Boniface, Mark)
# Cluster 3: 1, 2, 3, 6 (Romuald, Bonaventure, Ambrose, Louis)
# Cluster 4: 16, 17 (Basil, Elias)
# Notice that Cluster 2 is all "Young Turks", Cluster 4 is all 
#  "Outcasts", Cluster 3 is mostly "Loyal Opposition", and
#  Cluster 1 is a bit mixed. As BEJ note, this only used
#  a small portion of Sampson's data, but still did okay.

# Matrix 12.4
# Nothing to see here...

# Figure 12.3
# Not doing this one. Just note how the four "blocks" of
#  Matrix 12.4 become the matrix shown in Figure 12.3, and that
#  the network created from that "block matrix" comes about.
#  (The vertex on the left is block 1, while the vertex on the
#  right is block 2.)

# Figure 12.4
equiv.clust(esteem_mat) -> esteem_equiv
equiv.clust(disesteem_mat) -> disesteem_equiv

blockmodel(dat = esteem_mat, 
           ec = esteem_equiv,
           k = 4) -> esteem_bm   # k = 4 because BEJ did
# You can get the groups here; the matrix isn't printed out, though.
esteem_bm

blockmodel(dat = disesteem_mat, 
           ec = disesteem_equiv,
           k = 4) -> disesteem_bm   # k = 4 because BEJ did
disesteem_bm

# Can also do:
blockmodel(dat = sampson_arr,
           ec = equiv_res,
           k = 4) -> sampson_bm
sampson_bm
# None of these are exactly BEJ, but notice that there are many
#  approaches to creating the equivalnece, and I don't know which
#  one UCINET uses. Still, the concept is the same.

# Matrix 12.5
# The "number of changes" is also known as the "Hamming distance"
#  Hmm...by default, equiv.cluster() and sedist9) use the Hamming
#  distance.
blockmodel(dat = esteem_mat, 
           ec = esteem_equiv,
           k = 3) -> esteem_bm   # k = 3 because BEJ did
# You can get the groups here; the matrix isn't printed out, though.
esteem_bm
# Still not identical. See the comments for Matrix 12.4

# Table 12.1
# I don't think I have the data to do these, but it would use
#  those data if I did have them.

# Matrix 12.6
# No data, and it isn't that interesting to create the data.

# Matrix 12.7
fread(here("Data/TARO.csv"), header = TRUE) -> taro_mat
as.matrix(taro_mat[,-1]) -> taro_mat
colnames(taro_mat) -> rownames(taro_mat)
equiv.clust(taro_mat) -> taro_clust
blockmodel(taro_mat, taro_clust, k = 3) -> taro_bm
taro_bm
# Again, not identical, but that probably has to do with the 
#  clustering methods.

# Figure 12.5

# Figure 12.6
# I know what BEJ are trying to get at here, but this would be
#  a lot of work for a not-very-clear example.

# Matrix 12.8
fread(here("Data/Baker.csv"),
      header= TRUE) -> baker_mat
as.matrix(baker_mat[,-1]) -> baker_mat
colnames(baker_mat) -> rownames(baker_mat)
baker_mat

# Table 12.3
#Compute the k-cores based on total degree; the "discrete method"
kc<-kcores(baker_mat)
kc
# Sorry folks, I know nothing about "continuous core-periphery"
#  scores, so I can't do this. (I must say, however, that even
#  BEJ state that "[i]t follows that the continuous method, which
#  must take account of the off-diagonal blocks, cannot 
#  differentiate between core–periphery and periphery–core 
#  interaction. This means that it is better suited to symmetric
#  data (or at least data that are very nearly symmetric." So, 
#  probably not a big loss.)

# Figure 12.7
# Plot the result. gplot() plots a matrix as if it
#  were a graph.
gplot(baker_mat,
      vertex.col=kc)

## Chapter 13 - still to do ####

# Matrix 13.1
# Same as Matrix 2.3
load("Data/DeepSouth.RData")
as.matrix(davisDyn) -> davisDyn

# Matrix 13.2
davisDyn %*%  t(davisDyn) -> davis_women

# Figure 13.1
graph_from_adjacency_matrix(davis_women,
                            mode = "undirected",
                            weighted = TRUE,
                            add.colnames = NULL,
                            diag = FALSE) -> davis_women_gr
set.seed(42)
# This would plot the graph:
plot(davis_women_gr,
     shape = "square",
     edge.width = E(davis_women_gr)$weight)
# However, BEJ seem to have filtered out weights <= 2. Hence, we
#  induce the subgraph, and then plot:
subgraph.edges(davis_women_gr,
               eids = which(E(davis_women_gr)$weight > 2)) -> davis_2_gr
set.seed(42)
plot(davis_2_gr,
     vertex.shape = "square",
     edge.width = E(davis_2_gr)$weight)

# Matrix 13.3
t(davisDyn) %*% davisDyn -> davis_event
davis_event
graph_from_adjacency_matrix(davis_event,
                            mode = "undirected",
                            weighted = TRUE,
                            add.colnames = NULL) -> davis_event_gr

# Matrix 13.4
# Looking at the UCINET help (online) for "Data Affiliations" gives this:
# "Converts an m´n matrix to an m´m or n´n by forming AA'  or A'A using two 
#  different types of binary multiplication...."
# That's just what we did. However, they also say:
# 
# "The routine also allows for the final matrix to be normalized to accommodate
#  the different sizes of the events. Consider two actors i and j and let X be
#  the product of the number events they both attended and the number of events
#  they both did not attend, let Y be the product of the number events i 
#  attended and j did not with the number of events j attended and i did not. 
#  If X=Y the normalized entry is 0.5 otherwise it is (X-SQRT(XY))/(X-Y)."
#
# The last part is how to normalize this stuff. For the original matrix
#  Events both attended: sum(row(i) * row(j)) = both
#  Events both did not attend: number of columns - both
#  i attended without j = sum(row(i)) - both
#  j attended without i = sum(row(j)) - both
# After doing all this, then normalize to 100.

# Figure 13.2
# From Matrix 13.4

# Matrix 13.5
matrix(0, ncol = ncol(davisDyn), nrow = ncol(davisDyn)) -> lower_right
matrix(0, ncol = nrow(davisDyn), nrow = nrow(davisDyn)) -> upper_left
cbind(upper_left, davisDyn) -> upper
cbind(t(davisDyn), lower_right) -> lower
rbind(upper, lower) -> davis_bipart
davis_bipart

# Figure 13.3
graph_from_incidence_matrix(davisDyn,
                            add.names = NULL,
                            weighted = TRUE) -> davis_bipart_gr
set.seed(42)
plot(davis_bipart_gr,
     vertex.shape = ifelse(V(davis_bipart_gr)$type == TRUE,
                           "circle",
                           "square"))
# Figure 13.4
# I can't find enough detail to create a bi-clique anlaysis, so I can't
#  do the clustering.

# Figure 13.5
# install.packages("ITNr") - # Only one mode graphs, sadly
# library(ITNr)
# core_periphery_weighted(davis_bipart_gr,
#                         type = "undirected")

# Figure 13.6
# structural equivalence for women & for events

# Figure 13.7
# direct blockmodel for davisDyn
#

# Figure 13.8
# structural blockmodel for davisDyn

## Chapter 14 - still to do ####
# 5 Figs, 1 Mat, 2 Tables

# Figure 14.1
# Be careful...this takes a while to plot...so limit iterations
fread(here("Data/pv960.csv"), header = TRUE) -> pv960_mat

set.seed(42)
list("niter" = 10) -> layout_pars  # 10 iterations, not 500
# Doesn't look so good, but only 10 seconds
gplot(pv960_mat,
      gmode = "graph",
      vertex.sides = 4,
      vertex.rot = 45,
      layout.part = layout_pars)

# Figure 14.2
# Be careful...this takes a while to plot if you don't limit
#  the iterations for the Fruchterman-Reingold
as.matrix(fread(here("Data/pv960.csv"),
                header = TRUE)) -> pv960_mat

set.seed(42)
list("niter" = 10) -> layout_pars  # 10 iterations, not 500
# Doesn't look so good, but only a few seconds
gplot(pv960_mat,
      thresh = 2,           # Only ties >2 displayed
      gmode = "graph",
      vertex.sides = 4,
      vertex.rot = 45,
      displayisolates = FALSE,
      layout.par = layout_pars)

# Figure 14.3
# The figure caption says "after removing ties with
#  edge weights less than 3". So, let's do that first:
0 -> pv960_mat[which(pv960_mat<3)]

network(pv960_mat, directed = FALSE) -> pv_net
kcores(pv_net, mode = "graph") -> kc

# I could induce a subgraph of pv_net, but this is faster and involves
#  less typing.
which(kc > 10) -> kc10
pv960_mat[kc10, kc10] -> pv10
set.seed(42)
gplot(pv10,
      gmode = "graph")

# Table 14.1
# Your times may vary. There are several ways to "benchmark" the running
#  times of functions in R, like with package:microbenchmark.

# Matrix 14.1
# I'd love to do this, but without the department information, I
#  really can't.
read.ucinet.header(here("Data/pv504")) -> pv504_hdr
read.ucinet(here("Data/pv504")) -> pv504_mat

# Figure 14.4
# See the comment for Matrix 14.1

# Figure 14.5
# Sure, this is a network, but it is trying to help you understand
#  Louvain, so we won't reproduce it here.



## Chapter 15 ####

# Figure 15.1
# Data not available
# You can, however, get your own ego network from Facebook, using 
#  the package RFacebook, but as Facebook changes its access policies,
#  the package is sometimes in compliance and sometimes not.
#  See http://thinktostart.com/analyzing-facebook-with-r/ for more info
#  on collecting data.

# Table 15.1
# Nothing to see here...move along

# Table 15.2
# Nothing to see here...move along

# Figure 15.2
# Nothing to see here...move along

# Figure 15.3
# Nothing to see here...move along

# Table 15.3
# Nothing to see here...move along

# Table 15.4
# Nothing to see here...move along

# Matrix 15.1
# Nothing to see here...move along

# Figure 15.4
# Nothing to see here...by now you should know how to do this one.

# Table 15.5
# Nothing to see here...move along

# Figure 15.5
# Nothing to see here...move along

# Figure 15.6
# Nothing to see here...move along

# Figure 15.7
# Nothing to see here...move along

# Table 15.6
# Nothing to see here...move along

# Figure 15.8
# Nothing to see here...move along

# Table 15.7
# Nothing to see here...move along

# Figure 15.9
# Nothing to see here...move along

# Figure 15.10
# Nothing to see here...move along

# Figure 15.11
# Nothing to see here...move along

# Figure 15.12
# Nothing to see here...move along
