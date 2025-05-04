umap_coord_anno <- function(size) {
  ggplot(tibble(group = c("UMAP1", "UMAP2"),
                x = c(0, 0), xend = c(1, 0),
                y = c(0, 0), yend = c(0, 1),
                lx = c(0.5, -0.15), ly = c(-0.15, 0.5),
                angle = c(0, 90))) +
    geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
                 arrow = arrow(angle = 20, type = "closed", length = unit(0.1, "npc")),
                 size = size/4, lineend = "round") +
    geom_text(aes(lx, ly, label = group, angle = angle), size = size) +
    theme_void() +
    coord_fixed(xlim = c(-0.3, 1), ylim = c(-0.3, 1))
}

add_umap_coord <- function(gg_obj, x = 0, y = 0, width = 0.3, height = 0.3, size = 4.5) {
  p <- ggdraw() + 
    draw_plot(gg_obj, x = 0, y = 0, width = 1, height = 1) +
    draw_plot(umap_coord_anno(size = size), x = x, y = y, width = width, height = height)
  return(p)
}

remove_xaxis <- theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.line.x = element_blank())

remove_yaxis <- theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.title.y = element_blank(),
                      axis.line.y = element_blank())
