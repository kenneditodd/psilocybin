#' get_cluster_colors
#'
#' @param num_clusters
#'
#' @examples
get_cluster_colors <- function(num_clusters) {
  
  cluster_colors <- vector()
  
  if (num_clusters == 6) {
    cluster_colors <- c("#fa4e3e", "#fabf3e", "#30d150", "#3086d1", "#b186f7", 
                        "#fa78e9")
  } else if (num_clusters == 7) {
    cluster_colors <- c("#fa4e3e", "#fabf3e", "#30d150", "#78f3fa", "#2966f2", 
                        "#b186f7", "#fa78e9")
  } else if (num_clusters == 8) {
    cluster_colors <- c("#fa4e3e", "#fabf3e", "#55f531","#08780c", "#78f3fa", 
                        "#2966f2", "#b186f7", "#fa78e9")
  } else if (num_clusters == 9) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#55f531","#08780c", 
                        "#78f3fa", "#2966f2", "#b186f7", "#fa78e9")
  } else if (num_clusters == 10) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#55f531","#08780c", 
                        "#0ffcf8", "#117ea6", "#072be3", "#b186f7", "#fa78e9")
  } else if (num_clusters == 11) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#85fa3c","#3f940a", 
                        "#1b5212", "#0ffcf8", "#117ea6", "#072be3", "#b186f7", 
                        "#fa78e9")
  } else if (num_clusters == 12) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#85fa3c","#3f940a", 
                        "#1b5212", "#0ffcf8", "#117ea6", "#072be3", "#d996fa", 
                        "#ac11fa", "#ff75bc")
  } else if (num_clusters == 13) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#85fa3c","#3f940a", 
                        "#1b5212", "#0ffcf8", "#117ea6", "#072be3", "#d996fa", 
                        "#ac11fa", "#ff75bc", "gray")
  } else if (num_clusters == 14) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#85fa3c","#3f940a", 
                        "#1b5212", "#0ffcf8", "#117ea6", "#072be3", "#d996fa", 
                        "#ac11fa", "#ff75bc", "gray", "chocolate4")
  } else if (num_clusters == 15) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#85fa3c","#3f940a", 
                        "#1b5212", "#0ffcf8", "#117ea6", "#072be3", "#d996fa", 
                        "#ac11fa", "#ff75bc", "gray", "chocolate4", "black")
  } else if (num_clusters == 16) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#85fa3c","#3f940a", 
                        "#1b5212", "#0ffcf8", "#117ea6", "#072be3", "#d996fa", 
                        "#ac11fa", "#ff75bc", "gray", "tan", "chocolate4", 
                        "black")
  } else if (num_clusters == 17) {
    cluster_colors <- c("#fc4a0f", "#fca50f", "#faec28", "#85fa3c","#3f940a", 
                        "#1b5212", "#0ffcf8", "#117ea6", "#072be3", "#d996fa", 
                        "#a916f2", "#ffa1ef", "#fc1cd7", "gray", "tan", 
                        "chocolate4", "black")
  } 
  
  return(cluster_colors)
  
}
