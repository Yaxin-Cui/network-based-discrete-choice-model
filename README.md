# network-based-discrete-choice-model
Tutorial of network-based discrete choice modeling based on bipartite Exponential Random Graph Model (ERGM) approach. This demonstration contains code and data associated with the paper "A network-based discrete choice model for decision-based design."

## Notes

1. The data are anonymized, with each feature denoted by V1, V2, V3, and V4. The first three features (V1, V2, and V3) are continuous variables after scaling, while V4 is a categorical variable.
2. The code is divided into four parts:
   - Part 1: Load data and network construction
   - Part 2: Network visualization
   - Part 3: ERGM estimation
   - Part 4: ERGM prediction
3. The required package is the `statnet` package.

## Data

The data files should be placed in the same directory as this Markdown file or in a specified path. The data files are:

- `Sampled_data_to_share.csv`
- `test_data_to_share.csv`

## Part 1: Load Data and Network Construction

```R
# data1-5000 customer, each has 6 products in choice set and 1 product finally buy
data_train <- read.csv("Sampled_data_to_share.csv", header = TRUE) 
data_train$V4 <- as.factor(data_train$V4)
data_train$V4 <- relevel(data_train$V4, "A")

# function to set attributes
set_attr <- function(df, net, attr_product){
  df$src <- df$rspd_id
  df$dest <- df$model_id
  unique_rspd <- unique(df$rspd_id)
  unique_model <- unique(df$model_id)
  
  set.vertex.attribute(net, attrname = names(attr_product)[1:(ncol(attr_product)-1)], value = attr_product[1:(ncol(attr_product)-1)], v = attr_product$dest)
  return(net)
}

# function to create networks based on the data frame
# function to create networks based on the data frame
make_network <- function(df){
  # add src and dest id to dataframe
  df$src <- df$rspd_id
  df$dest <- df$model_id
  unique_rspd <- unique(df$rspd_id)
  unique_model <- unique(df$model_id)
  num_rspd <- length(unique_rspd)
  num_model <- length(unique_model)
  for(i in 1:nrow(df)){
    df$src[i] <- which(unique_rspd == df$src[i])
    df$dest[i] <- which(unique_model == df$dest[i]) + num_rspd
  }
  
  # add product attributes here
  attr_product <- unique(df[,c("V1", "V2", "V3", "V4", "dest")])
  
  ## inverse consideration network
  el = as.matrix(df[c('src','dest')])
  net_consideration <- network(x = el, matrix.type = "edgelist", directed = F, 
                               bipartite = num_rspd)
  mat_inv <- 1- as.matrix.network(net_consideration)
  
  
  ## purchase network
  el_purchase = as.matrix(df[df$purchase == 1, c('src','dest')])
  attr(el_purchase,'n') = length(unique_model)+length(unique_rspd)
  net_purchase <- network(x = el_purchase, matrix.type = "edgelist", directed = F,
                          bipartite = num_rspd)
  
  net_purchase <- set_attr(df, net_purchase, attr_product)
  net_consideration <- set_attr(df, net_consideration, attr_product)
  
  return(list("net_purchase" = net_purchase, "net_consideration" = net_consideration, "mat_inv" = mat_inv))
}

# make network lists: inclued net_purchase, net_consideration and mat_inv
newList <- make_network(data_train)

net_purchase <- newList$net_purchase
net_consideration <- newList$net_consideration
mat_inv <- newList$mat_inv

# check the summary statistics of net_purchase and net_consideration
summary(net_purchase)
summary(net_consideration)

