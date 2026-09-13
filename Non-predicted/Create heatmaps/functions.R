library(plotly)
library(webshot2)

#function to scale the number of cols of a matrix to a new length
scaler <- function(matrix, new_width){
  old_width <- ncol(matrix)
  #increase length (nearest neighbor resampling)
  if(old_width < new_width){
    #create vector of length new_width with equally spaced values from 1 to old_width
    idx <- round(seq(1, old_width, length.out=new_width))
    #create matrix with new_width
    return(matrix[,idx])
  } 
  #decrease length (block aggregation)
  #when old_width == new_width nothing changes here
  else{
    #divide the cols in #new_width blocks
    blocks <- cut(1:old_width, breaks=new_width, labels=FALSE)
    #go threw every block (cols) and calculate the max value per row (0 = no data, 1 = homozygous reference, 2 = heterozygous variant, 3 = homozygous variant)
    #Attention: NA values are possible (thats why I handle them in the last function)
    out <- sapply(split(1:old_width,blocks), function(cols){apply(matrix[,cols,drop=FALSE],MARGIN = 1, function(row){row <- row[!is.na(row)]; if(length(row)==0) return(NA) else max(row)})})
    return(out)
  }
}

#function to create a heatmap
create_heatmap <- function(matrix, base_one, miRNA_one, loop, miRNA_two, base_two, zmin, zmax, name, print_scale, colorscale, path){
  #define positions, colors and names for the different parts of the miRNA
  regions <- list(
    list(start=base_one, end=miRNA_one-1, color="rgba(200,200,200,0.2)", name="base"),
    list(start=miRNA_one, end=loop-1, color="rgba(150,150,160,0.2)", name="miRNA-5p"),
    list(start=loop, end=miRNA_two-1, color="rgba(120,120,120,0.2)", name="loop"),
    list(start=miRNA_two, end=base_two-1, color="rgba(150,150,160,0.2)", name="miRNA-3p"),
    list(start=base_two, end=base_two+24, color="rgba(200,200,200,0.2)", name="base"))
  
  region_shapes <- lapply(regions, function(r) {
    list(type = "rect", xref = "x", x0 = r$start - 0.5, x1 = r$end + 0.5, yref = "paper", y0 = 0, y1 = 1, fillcolor = r$color, line = list(width = 0), layer = "below")})
  
  region_annotations <- lapply(regions, function(r) {
    list(x = (r$start + r$end) / 2, y = 0, yref = "paper", text = r$name, showarrow = FALSE, yanchor = "top", yshift = -5, font = list(size = 9, color = "black"))})
  
  region_lines <- lapply(regions[2:5], function(r) {
    list(type = "line", xref = "x", x0 = r$start - 0.5, x1 = r$start - 0.5, yref = "paper", y0 = -0.01, y1 = 0, line = list(color = "black", width = 0.5))}) 
  
  #create heatmap
  p <- plot_ly(z = matrix, x = 1:(ncol(matrix)), type = "heatmap", colorscale = colorscale, zmin = zmin, zmax = zmax, showscale = print_scale) %>% 
    
    config(displayModeBar = FALSE) %>%
    
    layout(shapes = c(region_shapes, region_lines), annotations = region_annotations, xaxis = list(side = "top", showgrid = FALSE), yaxis = list(autorange = "reversed", showgrid = FALSE),
           legend = list(font = list(size = 8,family = "Arial")), title = list(text = paste0("<b>", name), y = 1, yanchor = "top", yref = "container", x = 0.52, font = list(size = 20,family = "Arial")), margin = list(t = 60))
  
  #save heatmap as html
  setwd(paste(path, "/html", sep = ""))
  htmlwidgets::saveWidget(as_widget(p), paste(name, ".html", sep = ""))
  #save heatmap as pdf
  webshot(url = paste(name, ".html", sep = ""), file =  paste(path, "/pdf/", name, ".pdf", sep = ""), delay = 2) 
}

#create heatmaps for the haplotypes (y = rownames(matrix) is needed)
create_heatmap_haplo <- function(matrix, base_one, miRNA_one, loop, miRNA_two, base_two, zmin, zmax, name, print_scale, colorscale, path){
  #define positions, colors and names for the different parts of the miRNA
  regions <- list(
    list(start=base_one, end=miRNA_one-1, color="rgba(200,200,200,0.2)", name="base"),
    list(start=miRNA_one, end=loop-1, color="rgba(150,150,160,0.2)", name="miRNA-5p"),
    list(start=loop, end=miRNA_two-1, color="rgba(120,120,120,0.2)", name="loop"),
    list(start=miRNA_two, end=base_two-1, color="rgba(150,150,160,0.2)", name="miRNA-3p"),
    list(start=base_two, end=base_two+24, color="rgba(200,200,200,0.2)", name="base"))
  
  region_shapes <- lapply(regions, function(r) {
    list(type = "rect", xref = "x", x0 = r$start - 0.5, x1 = r$end + 0.5, yref = "paper", y0 = 0, y1 = 1, fillcolor = r$color, line = list(width = 0), layer = "below")})
  
  region_annotations <- lapply(regions, function(r) {
    list(x = (r$start + r$end) / 2, y = 0, yref = "paper", text = r$name, showarrow = FALSE, yanchor = "top", yshift = -5, font = list(size = 9, color = "black"))})
  
  region_lines <- lapply(regions[2:5], function(r) {
    list(type = "line", xref = "x", x0 = r$start - 0.5, x1 = r$start - 0.5, yref = "paper", y0 = -0.01, y1 = 0, line = list(color = "black", width = 0.5))}) 
  
  #create heatmap
  p <- plot_ly(z = matrix, x = 1:(ncol(matrix)), y = rownames(matrix), type = "heatmap", colorscale = colorscale, zmin = zmin, zmax = zmax, showscale = print_scale) %>% 
    
    config(displayModeBar = FALSE) %>%
    
    layout(shapes = c(region_shapes, region_lines), annotations = region_annotations, xaxis = list(side = "top", showgrid = FALSE), yaxis = list(autorange = "reversed", showgrid = FALSE),
           legend = list(font = list(size = 8,family = "Arial")), title = list(text = paste0("<b>", name), y = 1, yanchor = "top", yref = "container", x = 0.52, font = list(size = 20,family = "Arial")), margin = list(t = 60))
  
  #save heatmap as html
  setwd(path)
  htmlwidgets::saveWidget(as_widget(p), paste0(name, ".html"))
  #save heatmap as pdf
  webshot(url = paste0(name, ".html"), file =  paste0(name, ".pdf"), delay = 2) 
}

