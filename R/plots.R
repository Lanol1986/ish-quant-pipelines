
ish_theme <- function(base_size=12) {
  ggplot2::theme_classic(base_size=base_size) +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill="white",colour=NA),
      plot.background=ggplot2::element_rect(fill="white",colour=NA),
      legend.position="bottom",strip.background=ggplot2::element_rect(fill="white",colour="grey40"),
      strip.text=ggplot2::element_text(face="bold"),axis.text=ggplot2::element_text(colour="black"),
      plot.title=ggplot2::element_text(face="bold"),
      plot.caption=ggplot2::element_text(colour="grey40",hjust=0))
}

ish_theme_presentation <- function(base_size=14) {
  ish_theme(base_size)+ggplot2::theme(plot.title=ggplot2::element_blank(),
    plot.subtitle=ggplot2::element_blank(),plot.caption=ggplot2::element_text(size=8))
}

ish_save_plot <- function(p,path,width=9,height=5.5,dpi=300) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  ggplot2::ggsave(path,p,width=width,height=height,dpi=dpi,bg="white")
  invisible(path)
}

ish_plot_panel_pair <- function(data,value_col,y_label,title=NULL,reference=NULL,presentation=FALSE) {
  p <- ggplot2::ggplot(data,ggplot2::aes(x=panel,y=.data[[value_col]],group=animal_id))+
    ggplot2::geom_line(colour="grey70",linewidth=.5)+
    ggplot2::geom_point(size=2.8,colour="#D4A017")+
    ggplot2::stat_summary(ggplot2::aes(group=panel),fun=mean,geom="crossbar",
                          width=.45,colour="black",linewidth=.45)+
    ggplot2::labs(x=NULL,y=y_label,title=title,
      caption="Panels are repeated consecutive-section measurements from the same animals and do not increase biological n.")
  if (!is.null(reference)) p <- p+ggplot2::geom_hline(yintercept=reference,linetype="dashed",colour="grey45")
  p + if (presentation) ish_theme_presentation() else ish_theme()
}

ish_plot_panel_facets <- function(data,value_col,x_col,y_label,title=NULL,reference=NULL,presentation=FALSE) {
  p <- ggplot2::ggplot(data,ggplot2::aes(x=.data[[x_col]],y=.data[[value_col]]))+
    ggplot2::geom_point(size=2.6,colour="#D4A017",
                        position=ggplot2::position_jitter(width=.06,height=0))+
    ggplot2::stat_summary(fun=mean,geom="crossbar",width=.4,colour="black",linewidth=.45)+
    ggplot2::facet_wrap(~panel)+
    ggplot2::labs(x=NULL,y=y_label,title=title,caption="One point = one animal; panel identity is preserved.")
  if (!is.null(reference)) p <- p+ggplot2::geom_hline(yintercept=reference,linetype="dashed",colour="grey45")
  p + if (presentation) ish_theme_presentation() else ish_theme()
}

ish_plot_cell_proportions <- function(comp) {
  ggplot2::ggplot(comp,ggplot2::aes(x=animal_id,y=pct_cells,fill=cell_type))+
    ggplot2::geom_col(width=.82)+
    ggplot2::facet_grid(group~panel,scales="free_x",space="free_x")+
    ggplot2::scale_y_continuous(limits=c(0,100),expand=ggplot2::expansion(mult=c(0,.02)))+
    ggplot2::labs(x=NULL,y="Cell proportion (%)",fill="Cell type",
      title="Cell composition by animal and panel",
      caption="ROI/segmentation composition QC; not an estimate of true whole-brain cell proportions.")+
    ish_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
}

ish_spatial_plot <- function(data,endpoint,cfg,title=NULL) {
  cols <- paste0("spatial__",c("xmin","xmax","ymin","ymax"))
  if (!all(cols %in% names(data))) stop("Spatial bounding-box columns missing.",call.=FALSE)
  d <- data |> dplyr::mutate(x_centroid=(spatial__xmin+spatial__xmax)/2,
                             y_centroid=(spatial__ymin+spatial__ymax)/2)
  lo <- cfg$spatial$percentile_lower %||% .01
  hi <- cfg$spatial$percentile_upper %||% .99
  lim <- as.numeric(stats::quantile(d[[endpoint]][is.finite(d[[endpoint]])],c(lo,hi),na.rm=TRUE,names=FALSE))
  d$display_value <- pmin(pmax(d[[endpoint]],lim[1]),lim[2])
  p <- ggplot2::ggplot(d,ggplot2::aes(x=x_centroid,y=y_centroid,colour=display_value))+
    ggplot2::geom_point(size=cfg$spatial$point_size %||% .35,alpha=cfg$spatial$alpha %||% .85)+
    ggplot2::scale_colour_viridis_c(option="magma",limits=lim,oob=scales::squish,name=endpoint)+
    ggplot2::scale_y_reverse()+ggplot2::coord_equal()+
    ggplot2::labs(x=NULL,y=NULL,title=title,
      caption=sprintf("Magma; %.0fth-%.0fth percentile display limits. Coordinates are cell centroids.",lo*100,hi*100))+
    ish_theme()+ggplot2::theme(axis.text=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank())
  list(plot=p,metadata=tibble::tibble(endpoint=endpoint,lower_percentile=lo,
    upper_percentile=hi,display_min=lim[1],display_max=lim[2],palette="magma"))
}

ish_presentation_mrna <- function(data,control_group) {
  d <- data |> dplyr::filter(endpoint_role=="target_mrna",metric=="mean_all",
    as.character(group)!=control_group,is.finite(normalized_value)) |>
    dplyr::mutate(display_value=pmin(pmax(normalized_value,0),100),
                  display_capped=normalized_value<0|normalized_value>100)
  ggplot2::ggplot(d,ggplot2::aes(x=cell_type,y=display_value))+
    ggplot2::geom_hline(yintercept=100,linetype="dashed",colour="grey45")+
    ggplot2::geom_point(data=dplyr::filter(d,!display_capped),size=2.7,colour="#D4A017")+
    ggplot2::geom_point(data=dplyr::filter(d,display_capped),shape=17,size=3,colour="#D4A017")+
    ggplot2::stat_summary(fun=mean,geom="crossbar",width=.4,colour="black",linewidth=.45)+
    ggplot2::facet_wrap(~panel)+ggplot2::coord_cartesian(ylim=c(0,100))+
    ggplot2::labs(x=NULL,y="Target mRNA remaining (% vehicle)",
      caption="Triangle = source value outside the 0-100% display range; stored values remain uncapped.")+
    ish_theme_presentation()
}
