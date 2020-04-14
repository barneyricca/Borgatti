# The following script is a "port" of Borgatti, Everett, & Johnson (2nd
#  edition) to R.
# In a couple locations, I have used fread() to open ".dat" files, which is
#  a common network file format. I have also included a couple functions to
#  read basic UCINET files, which are also a common format. However, I think
#  that all the data are also included as .CSV files; I ported the data from
#  various sources.
#
# If you installed the package from GitHub, you should have the data as well.

## Preparation ####

#
# There are several other potentially useful packages, including:
#   keyplayer - advanced centality methods
#   multiplex - for multplex graphs
#   ndtv - Network Dynamic Temporal Visualization
#   RSiena - particularly useful for longitudinal networks; see Snijders' 
#            Siena home page
#   sand - from a more mathematical text
#   tnet - miscellaneous advanced analyses
#   vegan - Some useful network and clustering functions
# They are mostly beyond the scope of this course, however, and so are not
#  included here.

{
  c("ca",            # Correspondence analysis
    "conflicted",    # To deal with conflicting function names
                     # I've had some strangeness with this
                     #  script. I suspect package:conflicted,
                     #  but I don't yet know for sure.
    "data.table",    # Fast data input/output
    "dplyr",         # This is all of tidyverse that gets used here
    "dtplyr",        # dplyr syntax with a data.table backend
#  "ggnetwork",      # ggplot-like plotting for networks
#  "graphlayouts",   # For plotting stuff
    "here",          # To find/store files w/o setwd() and getwd()
    "igraph",        # Basic network tools; we'll use statnet mostly
    "igraphdata",    # Some useful datasets
    "intergraph",    # Translate between igraph and statnet formats
#  "qgraph",         #
    "statnet"        # A suite of network tools, including ERGM and more
#  "tidyverse"       # Data manipulation
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

## A Note about plotting ####
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
netset.edge.value(monks_net, "Weight", monks_mat)

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

# Figure 6.4
length(network.vertex.names(F6.3_net)) -> num_vert
matrix(0,
       nrow = num_vert,
       ncol = 4,
       dimnames = list(network.vertex.names(F6.3_net),
                       c("Degree", "Closeness",
                         "Betweenness", "Eigenvector"))) -> cent_mat
# degree - must scale to get BEJ Table 6.2
100 * sna::degree(F6.3_net, gmode = "graph") / (num_vert-1) -> cent_mat[,"Degree"]

# closeness - must convert to percent to get BEJ Table 6.2
# I get different values no matter how I rescale; I think that BEJ are wrong.
#  Or maybe I'm wrong.
100 * sna::closeness(F6.3_net, gmode = "graph") -> cent_mat[,"Closeness"]

# betweenness - must scale to get BEJ Table 6.2
# I get different values no matter how I rescale; I think that BEJ are wrong.
#  Or maybe I'm wrong.
sna::betweenness(F6.3_net, gmode = "graph") / num_vert -> cent_mat[,"Betweenness"]

# eigenvector - must scale to get BEJ Table 6.2
# I get different values no matter how I rescale; I think that BEJ are wrong.
#  Or maybe I'm wrong.
sna::evcent(F6.3_net, gmode = "graph") -> cent_mat[,"Eigenvector"]

# Look at package:vegan because the matrix is not symmetric; correspondance
#  analysis
# Plot the points
# Plot the centrality scores too

# Figure 6.4

# Matrix 6.2 is simply city_mat above.
# Matrices 6.3, 6.4, 6.5 and 6.6 won't be redone here. However, here is the
#  clustering output of Figure 6.5:

hclust(city_dist, method="single") -> hc_single  # The BEJ way
plot(hc_single)
hc_single$height    # Read the dendogram from the bottom to see these in
                    #  context.

hc_singe$labels
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
#  investigate that problem right now.
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
camp_net %s% which(network::get.vertex.attribute(camp_net, "Sex") == 1) -> camp_red
set.seed(42)
plot.network(camp_red,
             label = network.vertex.names(camp_red),
             vertex.cex = sqrt(network::get.vertex.attribute(camp_net, "Betweenness")),
             pad = 1)

# I'm lazy...I'll skip 7.8 and 7.9

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
             label = network.vertex.names(MrHi_net))

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

layout_with_stress(davis_gr,
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

layout_with_stress(davis_red,
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

layout_with_stress(davis_gr,
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

# Figure 8.2


# From:
#  https://github.com/kateto/Network_Analysis_R_Examples/blob/master/R%20Scripts/Comm645-MRQAP.R
# Multiple regression for network variables: MRQAP
#=============================================================#


# Krackhardt data


# We can similarly use permutaton tests for a linear regression
# on network variables - known as MRQAP (multiple regression QAP)

# Note that the matrices used in the regression can contain different relations, the same
# relation at different points in time, or even dyadic covariate (attribute) data. For instance, 
# you can have a matrix of social ties as a DV and a same-gender attrubute matrix as an IV.
# (i.e. a matrix where (i,j)=1 if i and j are both female(male), and 0 if they have different genders)

# Let's create an array containing 3 random networks, each with 10 nodes:

x <- rgraph(10, 3) 
x[1,,] # matrix 1
x[2,,] # matrix 2
x[3,,] # matrix 3

# A matrix constructed as a linear combination of the first two networks in x:

y <- 3*x[1,,] + 5*x[2,,] + 7

# Now let's run a network regression of y on x and see what happens. 
# netlm takes a single matrix or network object as its first argument (DV) 
# and a stack of matrices (array containing all the IVs) or a list of 
# network objects as its second argument.

net.mod.1 <- netlm(y, x, reps=100)
summary(net.mod.1)

# Oh look - the first two parameters are significant and close to 3 and 5, and the intercept is 7
# Just as you'd expect based on the way we constructed y.


# Figure 8.3

# Load Newcomb's Fraternity Data
array(0,
      dim = list(15, 17, 17)) -> newcomb_arr
for(i in 1:15) {
as.matrix(fread("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/newfrat.dat",
                skip = (3 + i * 17),
                nrows = 17
                )) -> new_array[i,,]
}
# If the data aren't there, can get them this way:
#  load(file = "Data/newcomb.RData")

set.seed(24322)             # This is different than a UCINET seed of 24322

qaptest(newcomb_arr,
        FUN = gcor,
        g1 = 1,
        g2 = 2,
        reps = 50000) -> newcomb_qap   # BEJ do 50k reps

# Figure 8.4

# Campnet

# Matrix 8.1


## Chapter 9 ####


## Chapter 10 ####


## Chapter 11 ####


## Chapter 12 ####


## Chapter 13 ####


## Chapter 14 ####





# Lazega Lawyers Data Set
# This is part of the data set that Adhikari & Dabbs (2018) use. (Only the 36 partners are included; the 35 associates are not. There are, apparently, some ongoing privacy issues with the data set.) First, let's look at the data set:
data(lazega)               # Data in the format of the older igraph package.
if(graph_version(lazega) != "0.8.0") {
  upgrade_graph(lazega) -> lazega      # Upgrades to version 0.8.0
}
head(elist.lazega)
str(v.attr.lazega)

# Look at the help on lazega to see how the data are coded.

# Now, these data are already a network, so we don't yet have to create a network. Let's explore some plotting:

plot(lazega)     # Not very helpful
# Following A & D
layout_with_kk(lazega) -> lazega_kk_layout
plot(lazega, layout = lazega_kk_layout)    # Laid out, but not colored
{
  plot(lazega, layout = lazega_kk_layout,    # Layout and colored by status
       vertex.color = vertex_attr(lazega)$Status,
       vertex.label = NA)
  legend("topright", legend = c("Partner", "Associate"), 
         fill = categorical_pal(2))
}

{                                            # But I like (better) labels
  plot(lazega, layout = lazega_kk_layout,    # Layout and colored by status
       vertex.color = vertex_attr(lazega)$Status,
       vertex.size = 20,
       vertex.label.cex = 0.7)
  legend("topright", legend = c("Partner", "Associate"), 
         fill = categorical_pal(2))
}

# Do college alumni link together?
{                                            # But I like (better) labels
  plot(lazega, layout = lazega_kk_layout,    # Layout and colored by status
       vertex.color = vertex_attr(lazega)$School,
       vertex.size = 20,
       vertex.label.cex = 0.7)
  legend("topleft", legend = c("Harvard or Yale", "UConn", "Other"), 
         fill = categorical_pal(3))
}

# Much more can be done, as well; this is just some playing around with plots.

# Let's put some vertex data in.
set_vertex_attr(camp_gr, 
                "Sex", 
                # index=V(camp_gr), 
                value = camp_sex$Sexo) -> camp_gr
# And use this in a plot
c("circle", "square") -> shapes
plot(camp_gr, 
     vertex.shape = shapes[V(camp_gr)$Sex])

# Now, for the statnet version
network(camp_mat) -> camp_net  # statnet/sna version
gplot(camp_net)

network::set.vertex.attribute(camp_net, "Sex", camp_sex$Sexo)
# camp_sex$role -> camp_net %v% "Role"
# The previous line is no longer equivalent to:
network::set.vertex.attribute(camp_net, "Role", camp_sex$role)

gplot(camp_net,
      vertex.cex = 2.5,
      vertex.col = camp_net %v% "Role",
      vertex.sides = (camp_net %v% "Sex")+2)


# Or, using intergraph
asNetwork(camp_gr) -> camp2_net  # Transfers vertex and edge attributes
plot(camp2_net)
gplot(camp2_net,
      vertex.cex = 2.5,
      vertex.col = camp2_net %v% "Role",
      vertex.sides = (camp2_net %v% "Sex")+2)


# And, if you are really serious about needing to plot, you should either use
# [ggnetwork](https://cran.r-project.org/web/packages/ggnetwork/vignettes/ggnetwork.html)
# or tkplot, which allows you to manual adjust things. (But, it is rather buggy and twitchy, so you need to really want to do this before I would recommend it.)



# Some made-up data for F and E
matrix(c(0, 0, 0, 1, 0,
         0, 0, 0, 1, 0,
         1, 0, 0, 0, 1,
         1, 1, 0, 0, 0,
         0, 0, 1, 0, 0),
       nrow = 5,
       ncol = 5,
       byrow = TRUE) -> Fr

matrix(c(0, 0, 0, 0, 1,
         0, 0, 0, 0, 0,
         0, 1, 0, 0, 0,
         0, 0, 0, 0, 1,
         0, 0, 0, 1, 0),
       nrow = 5,
       ncol = 5,
       byrow = TRUE) -> En

Fr %*% En
Fr * En        # Incorrect!!!!!
En %*% Fr
Fr %*% Fr
Fr - En -> Re
Re


# Let's figure out how these work, exactly. (See pp. 23-25 of BEJ). _Fr %*% Fr_ is the number of times that _i_ has a friend that has _j_ as a friend. We can (and will, in BEJ 8,) determine if these things are more likely than we would expect by chance; in other words, we can do some inference.

# Harrison White (and others) looked at marriage rules and descent rules in various societies. In one scenario, he comes up with a number of interesting ideas and network analyses, encapsulated in 8 axioms:
  
1. _n_ clans in society
2. Permanent marriage rule (what clan marries what other clan)
3. Men from two different clans cannot marry women of the same clan
4. All children are assigned to a single clan (determined by mother and father)
5. Children whose fathers are in different clans must themselves be in different clans
6. A man can never marry a woman of his own clan.
7. Every person has at least one relative in every other clan.
8. "Whether two people who are related by marriage and descent links are in the same clan depends only on the kind of relationship, not on the clan either one belongs to".
(White, as quoted in Bradley & Meek, 1986.)

For example, Bradley & Meek use 5 clans (P, Q, R, S, and T), and the matrices:
  \[
    M = marriage = \begin{equation}
    \begin{bmatrix}
    0 & 1 & 0 & 0 & 0 \\ 
    0 & 0 & 1 & 0 & 0 \\ 
    1 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 1 \\
    0 & 0 & 0 & 1 & 0 \end{bmatrix} 
    \end{equation}
    \]

\[
  D = descent = \begin{equation}
  \begin{bmatrix}
  0 & 0 & 0 & 0 & 1 \\
  1 & 0 & 0 & 0 & 0 \\
  0 & 1 & 0 & 0 & 0 \\ 
  0 & 0 & 1 & 0 & 0 \\ 
  0 & 0 & 0 & 1 & 0 \end{bmatrix} 
  \end{equation}
  \]
In the marriage matrix, the rows (columns) are the husband's (wife's) clan; in descent matrix, the rows (columns) are the father's (child's) clan.

So, we can use this to find answer questions like what is the clan of my mother's brother's daughter? We do this by:
  
  * Transpose (i.e., exchange the rows and columns of) the descent matrix to find my father's clan. For matrices like we have, known as permutation matrices, the inverse and the transpose are the same operation. Inverting a matrix can be tricky, but we have R, so _solve(D)_ gives the inverse of C. _t(D)_ gives the transpose.
* Use the marriage matrix to find my mother's clan
* That clan is my mother's brother's clan
* Use the descent matrix to find my mother's brother's daughter's clan

In matrix notation, $t(D) %*% M %*% D$ gives the answer. Let's do that.
```{r MBD}
matrix(c(0, 1, 0, 0, 0,
         0, 0, 1, 0, 0,
         1, 0, 0, 0, 0,
         0, 0, 0, 0, 1,
         0, 0, 0, 1, 0),
       nrow = 5,
       ncol = 5,
       byrow = TRUE) -> M

matrix(c(0, 0, 0, 0, 1,
         1, 0, 0, 0, 0,
         0, 1, 0, 0, 0,
         0, 0, 1, 0, 0,
         0, 0, 0, 1, 0),
       nrow = 5,
       ncol = 5,
       byrow = TRUE) -> D

t(D) %*% M %*% D

```
The row is my clan, while the column is my MBD's clan.

Sometimes, we see some interesting block structure to these matrices. (And more often, we see some approximately block structure in the matrices; these approximate blocks are called communities.)

We can also have strength of ties, using numbers larger than one. For example, we could use for entries the number of times _i_ called _j_ in a week; those numbers would all be non-negative, but may be larger than 1.

Because of their usefulness, almost all computations on networks are actually done with matrices.

__Class 4: 23 January 2019__
Two "quick" review questions:

1. In the following network, identify one path, one trail that is NOT a path, and one walk that is neither a path nor a trail:
set.seed(42)
erdos.renyi.game(10, .3, directed = TRUE) -> ex_gr
plot(ex_gr)

2. For the example marriage and descent networks from the previous class, if you are in the Q clan, in which clan is your MDBS?
D

_Clustering_ is important, as it attempts to walk a middle ground between group-level properties (too course-grained) and individual-level properties (too fine-grained). We "[assign] items into groups or classes based on similarities or distances between them" (BEJ $\S$ 6.4). And while we do this all the time -  schools put students in tracks, and we all categorize people and events - complexity reduction can be problematic:

* Arrogance of reducer
* Violence done to the reduced

Still, let's do some of this. First, you should know that there are several uses of the word "cluster", so a modifier is always appropriate. BEJ 6 talks of _hierarchical clustering_, and gives the process; their example of distances between US cities is instructive (in Figure 6.5 remember the "Level" reads down the column; across that row is just the number that lists the oder in which the cities were originally given). We don't care so much about the process pseudocode as we do about getting an answer, so, let's repeat that example in R:

  # As BEJ
  # Plot the MDS!
  

Compare the two methods.

We can do this for most any data set.

head(mtcars)
dist(mtcars) -> mtcars_d
hclust(mtcars_d, method="average") -> cars_hc
plot(cars_hc)
str(cars_hc)


OK, let's look carefully at the output of hclust.

What else can we cluster? Here's some to look at:
  [Cluster examples](https://cran.r-project.org/web/packages/dendextend/vignettes/Cluster_Analysis.html)

Supporting Teachers
Once upon a time, I was on a project that surveyed beginning and experienced teachers to see how they were different. 

The theory (at the time) said that there were 9 things that teachers worried about: students, parents, colleagues, supervisors, preparation, autonomy, appearance, grading and workload. We, however, found only 8 factors: workload, appearance, colleagues, parents, grading, preparation, autonomy, and a "classroom" factor.

What are factors? Well, they are "a small number of *linear combinations* of the variables so as to capture most of the variation...With a large number of variables it may be easier to consider a small number of combinations of the original data rather than the entire dataframe." (Crawley, Chapter 25.) This relies on the data being somewhat correlated.

So, a 72-question survey was developed, and administered. The results are in "Supporting Teacher Data.csv" in the usual place. First, let's do a factor analysis on those data to see what the factors are. (The column headings are the theoretical category followed by the actual question number. E.g., Q0544 would be the 44th question, which theoretically applied to category 5.)

The experience question, although coded as numbers, really mean this:
Pre-service teachers:
1 - no methods courses & no observations
2 - methods courses but no observations
3 - observations but no student teaching
4 - currently student teaching
5 - other

Post:
1 - 1 to 5 years
2 - 6 to 10 years
3 - 11 to 15 years
4 - 16 to 20 years
5 - 21 or more years


library(data.table)
fread("http://citadel.sjfc.edu/faculty/bricca/Data/Beginning Teacher Factors(Z).csv")->fac.df

# Add rownumbers

fac.df %>%
  filter(.,group=="pre-service") %>%
  select(.,factor_classroom, factor_workload, factor_appearance, factor_colleagues,
         factor_parents, factor_grading, factor_preparation, factor_autonomy) -> temp.df
paste("T",1:nrow(temp.df),sep="")->labels
hclust(dist(temp.df))->hteach # Do this clustering first
plot(hteach,labels=labels, main= "")

From here, we see if we can identify what is important about the clusters, often by doing significance testing between the clusters.


How many clusters to use?
Some answers:
a) fvis() in the factoextra package
b) topology: Essentially, long lines on the dendograms
c) [Others that I'm not familiar enough with](http://www.rpubs.com/s_ritesh/Deciding_Clusters) to comment on except to say that they sometimes differ.
d) [Still others](https://www.statmethods.net/advstats/cluster.html)

What else can we cluster?
  
  HW 1: See Blackboard

__Class 5: 28 January 2019__
More on [How many clusters?](https://towardsdatascience.com/10-tips-for-choosing-the-optimal-number-of-clusters-277e93d72d92)

```{r clusters}
c(
  "magrittr",
  "cluster",
  "cluster.datasets",
  "cowplot",
  "NbClust",
  "clValid",
  "ggfortify",
  "clustree",
  "dendextend",
  "factoextra",
  "FactoMineR",
  "corrplot",
  "GGally",
  "ggiraphExtra",
  
  
) -> package_names  

data("all.mammals.milk.1956")
#raw_mammals <- all.mammals.milk.1956
# subset dataset
all.mammals.milk.1956 %>%    
  select(., -name) %>%       # subset data
  as_tibble(.) -> mammals    # Make it a tibble
# mammals <- as_tibble(mammals)

summary(mammals) %>% kable() %>% kable_styling()

mammals %>% 
  gather(Attributes, value, 1:5) %>% 
  ggplot(aes(x=value)) +
  geom_histogram(fill = "lightblue2", color = "black") + 
  facet_wrap(~Attributes, scales = "free_x") +
  labs(x = "Value", y = "Frequency")

corrplot(cor(mammals), type = "upper", method = "ellipse", tl.cex = 0.9)

mammals_scaled <- scale(mammals)
rownames(mammals_scaled) <- all.mammals.milk.1956$name

res.pca <- PCA(mammals_scaled,  graph = FALSE)
# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

# Extract the results for variables
var <- get_pca_var(res.pca)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# Control variable colors using their contributions to the principle axis
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) + theme_minimal() + ggtitle("Variables - PCA")

kmean_calc <- function(df, ...){
  kmeans(df, scaled = ..., nstart = 30)
}
km2 <- kmean_calc(mammals_scaled, 2)
km3 <- kmean_calc(mammals_scaled, 3)
km4 <- kmeans(mammals_scaled, 4)
km5 <- kmeans(mammals_scaled, 5)
km6 <- kmeans(mammals_scaled, 6)
km7 <- kmeans(mammals_scaled, 7)
km8 <- kmeans(mammals_scaled, 8)
km9 <- kmeans(mammals_scaled, 9)
km10 <- kmeans(mammals_scaled, 10)
km11 <- kmeans(mammals_scaled, 11)
p1 <- fviz_cluster(km2, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 2") 
p2 <- fviz_cluster(km3, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 3")
p3 <- fviz_cluster(km4, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 4")
p4 <- fviz_cluster(km5, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 5")
p5 <- fviz_cluster(km6, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 6")
p6 <- fviz_cluster(km7, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 7")
plot_grid(p1, p2, p3, p4, p5, p6, labels = c("k2", "k3", "k4", "k5", "k6", "k7"))



gap_stat <- clusGap(mammals_scaled, FUN = kmeans, nstart = 30, K.max = 24, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")



__Class 6: 30 January 2019__
Network stuff:
  
data(florentine)
plot(flobusiness)
plot(flomarriage)

array(dim=c(2,16,16)) -> florence
as.matrix.network(flobusiness,
                  matrix.type = "adjacency") -> florence[1,,]
as.matrix.network(flomarriage,
                  matrix.type = "adjacency") -> florence[2,,]
gcor(florence, g1=1, g2=2)


netlm(flobusiness, flomarriage) -> nl
summary(nl)

netlm(flobusiness, flomarriage, intercept = FALSE) -> nl1
summary(nl1)





#Generate three graphs
set.seed(20190130)
g<-array(dim=c(3,10,10))
g[1,,]<-rgraph(10)
g[2,,]<-rgraph(10,tprob=g[1,,]*0.8)
g[3,,]<-1; g[3,1,2]<-0              #This is nearly a clique

as.network(g[1,,]) -> n1
as.network(g[2,,]) -> n2
as.network(g[3,,]) -> n3

plot(n1)
plot(n2)
plot(n3)



#Perform qap tests of graph correlation
q.12<-qaptest(g,gcor,g1=1,g2=2, reps=1e3)
q.13<-qaptest(g,gcor,g1=1,g2=3)

#Examine the results
summary(q.12)
plot(q.12)
summary(q.13)
plot(q.13)



# qap tests of linear models
netlm(flobusiness, flomarriage, nullhyp = "qap")




Hypothesis Testing (Chapter 8); more ERGM?
  Permutation hypotheses
More stuff (dyadic hypotheses, node-level hypotheses, whole-network hypotheses)


Hypothesis Testing On Networks (BEJ 8)
__Class 7: 4 February 2019__
# ERGM
BEJ 8: What's important about ERGM? $\S$8.7, first paragraph.

Let's look at some examples

[ERGM Tutorial](https://michaellevy.name/blog/ERGM-tutorial/)

[ERGM Vignettes](https://cran.r-project.org/web/packages/ergm/vignettes/ergm.pdf)

[Another tutorial](https://www.win.tue.nl/~rmcastro/tmp/YES_VI/files/ERGM.pdf) - pp 23-27 for what is new/helpful here







Node characteristics
Edge characteristics
Layout (automatic, and adjusted; add cluster characteristics when clusters are done)


Cohesion
Reciprocity
Transitivity and clustering
Triads
Centralization & core-periphery indices



keyplayer package
Types of centrality


project topics due
cliques
algorithms
factions
modularity
directed networks, 



structural
blockmodels



Data Collection & Management (Chapters 4 & 5); contagion? API and scraping (e.g., scrape the Fisher Athletics rosters for name, sport, and major(s))




For matrix 5.7:
  # For the first, we need to know how many authors, and how
  #  many links are in each department.
  rowSums(scientists960_mat) - 
  diag(scientists960_mat) -> num_co_authors

data.frame("Author" = attributes960_df$V1,
           "Links" = num_co_authors,
           "Department" = attributes960_df$DeptID) -> authors_df

table(authors_df$Department) -> dept_numbers
authors_df %>%
  group_by(., Department) %>%
  summarise(.,
            "Papers per Department" = sum(Links)) -> dept_authors
c(1, 2, 3, 4, 5, 6, 7, 10, 20, 30) -> dept_ID
list() -> dept_members
for(i in 1:length(dept_ID)) {
  which(authors_df$Department == dept_ID[i]) -> dept_members[[i]]
}

dept_numbers
dept_authors
dept_members

rep(0,10) -> sum1
for(i in 1:10) {
  sum(scientists960_mat[1,dept_members[[i]]]) -> sum1[i]
}


# Because there are 10 departments, the latter for each author
#  is just the following:

# The number of times a department is a co-author is this:
(colSums(scientists960_mat) - 
    diag(scientists960_mat)) / 10 -> co_authored_per_dept




authors_df %>%
  group_by(., Department) %>%
  summarise(.,
            "Total Co-authors" = sum(LinksOut),
            "Total Co-authored" = sum(LinksIn)) -> links
