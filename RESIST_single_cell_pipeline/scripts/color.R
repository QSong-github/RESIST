library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library(ggsci)
library(R.utils)
library(scales)
library(tibble)

col<- c( "#E377C2FF", 
          '#FF7F0EFF',
         '#2CA02CFF', 
         "#9467BDFF",
         "#8C564BFF",
         "#FF9896FF",
        "#BCBD22FF" ,
        "#17BECFFF",
        "#AEC7E8FF",
        "#FFBB78FF",
        "#98DF8AFF"  ,
        "#C5B0D5FF" ,
        "#C49C94FF",
        '#F7B6D2FF',
        '#DBDB8DFF',
        '#9EDAE5FF') ## annotation color





col1<- c( '#E4C66F', '#5B9BD5' ) ## condition color

