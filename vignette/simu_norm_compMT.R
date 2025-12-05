#install.packages("aplpack")

library(MASS)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(png)
library(grid)
set.seed(1234)

save_dir = './figures/'

p=0.05
n=5000
theta=rbinom(n,1,p)
x=(1-theta)*rnorm(n,100)+theta*rnorm(n,103)
y=(1-theta)*rnorm(n,300)+theta*rnorm(n,298)
dat_normal=cbind(x,y)

x=rchisq(n,10)
y=rchisq(n,10)
dat_chisq=cbind(x,y)

# t-distribution settings
x = rt(n, df = 10)
y = rt(n, df = 10)
dat_t = cbind(x, y)

# log normal distribution settings
x = rlnorm(n, meanlog = 0, sdlog = 0.5)
y = rlnorm(n, meanlog = 0, sdlog = 0.5)
dat_lognormal = cbind(x, y)

# cauchy distribution settings
x = rcauchy(n, location = 0, scale = 1)
y = rcauchy(n, location = 0, scale = 1)
dat_cauchy = cbind(x, y)

dist_settings = c('compMT')
# dist_settings = c('chisq')
# dist_settings = c('t_dist')
# dist_settings = c('log_normal')
# dist_settings = c('cauchy')
inlier_settings = c(FALSE)
outter_settings = c(TRUE)
redefine_loop_settings = c(TRUE)
conservative_lambda_settings = c(1)
center_type_settings = c('hdepth')

all_settings = expand.grid(dist=dist_settings, normal_inlier=inlier_settings, normal_outter=outter_settings, conservative_lambda=conservative_lambda_settings, redefine_loop=redefine_loop_settings, center_type=center_type_settings)

for (setting_idx in 1:nrow(all_settings)) {
    dist_type = all_settings[setting_idx, 'dist']
    normal_inlier = all_settings[setting_idx, 'normal_inlier']
    normal_outter = all_settings[setting_idx, 'normal_outter']
    conservative_lambda = all_settings[setting_idx, 'conservative_lambda']
        redefine_loop = all_settings[setting_idx, 'redefine_loop']
        center_type = all_settings[setting_idx, 'center_type']
        fig_title = paste0(dist_type,
                                ifelse(center_type == 'hdepth', '_center_hdepth', '_center_mcd')
                                )
    dat = switch(as.character(dist_type),
                  'compMT' = dat_normal,
                  'chisq' = dat_chisq,
                  't_dist' = dat_t,
                  'log_normal' = dat_lognormal,
                  'cauchy' = dat_cauchy)
    bgplts = list()

    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='FWER',q=0.05, normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
            show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)

    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='FWER',q=0.1,normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
            show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)


    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='FWER',q=0.2, normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
            show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)

    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='FDR',q=0.01, normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
            show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)

    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='FDR',q=0.05,normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
            show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)


    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='FDR',q=0.1, normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
            show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)


    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='PFER',q=0.5, normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
                show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)

    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='PFER',q=1,normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
            show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)


    old=Sys.time()
    tmp_plt = bagplot(dat,factor=3,type1='PFER',q=5, normal_inlier=normal_inlier,normal_outter=normal_outter,conservative_lambda=conservative_lambda,create.plot=TRUE,approx.limit=10000, redefine_loop=redefine_loop, cex = 0.6, pch=1,
            show.outlier=TRUE,show.loophull=FALSE,
            show.bagpoints=TRUE,dkmethod=1,center_type = center_type,
            show.whiskers=TRUE,asymp_dist_pv = 'chisq', whisker.end.prop = 0.7,
            show.baghull=TRUE,verbose=FALSE,col.baghull = "#D3D3D3",
            show.fence_mag_bag = TRUE)
    bgplts = c(bgplts, list(tmp_plt))
    print(Sys.time()-old)

    plts_of_bgplts = lapply(bgplts, function(item) item$bgplt)

    # compute union xlim/ylim across all plots
    get_xy_limits <- function(p) {
      gb <- ggplot_build(p)
      xs <- unlist(lapply(gb$data, function(d) d[["x"]]), use.names = FALSE)
      ys <- unlist(lapply(gb$data, function(d) d[["y"]]), use.names = FALSE)
      xs <- xs[is.finite(xs)]
      ys <- ys[is.finite(ys)]
      list(
        x = if (length(xs)) range(xs, na.rm = TRUE) else c(NA_real_, NA_real_),
        y = if (length(ys)) range(ys, na.rm = TRUE) else c(NA_real_, NA_real_)
      )
    }
    lims_list <- lapply(plts_of_bgplts, get_xy_limits)
    xlim_union <- range(unlist(lapply(lims_list, function(l) l$x)), na.rm = TRUE)
    ylim_union <- range(unlist(lapply(lims_list, function(l) l$y)), na.rm = TRUE)

    # update each subplot with the union limits
    plts_of_bgplts <- lapply(plts_of_bgplts, function(p) {
      p + coord_cartesian(xlim = xlim_union, ylim = ylim_union) +
      theme(plot.margin = margin(t = 3.5, r = 1.3, b = 1.8, l = 1.8, unit = "lines"))
    })

    for (j in seq_along(plts_of_bgplts)) {
      ggsave(
        filename = file.path(
          save_dir,
          paste0("simu_", dist_type, "/", fig_title, "_sub", j, ".png")
        ),
        plot = plts_of_bgplts[[j]],
        width = 4,
        height = 4,
        bg = "white"
      )
    }
    
    subdir <- file.path(save_dir, paste0("simu_", dist_type))
    # apl_path  <- file.path(subdir, paste0(fig_title, "_aplpack.png"))
    sub_paths <- file.path(subdir, paste0(fig_title, "_sub", 1:9, ".png"))
    
    # Read images
    read_img_grob <- function(path) {
      img <- readPNG(path)
      rasterGrob(img, width = unit(1,"npc"), height = unit(1,"npc"))
    }
    # grobs <- c(list(read_img_grob(apl_path)), lapply(sub_paths, read_img_grob))
    grobs <- lapply(sub_paths, read_img_grob)
    
    # Wrap with labels
    labels <- c(
      "(a) FWER, q=0.05", "(b) FWER, q=0.1", "(c) FWER, q=0.2",
      "(d) FDR, q=0.01", "(e) FDR, q=0.05", "(f) FDR, q=0.1",
      "(g) PFER, q=0.5", "(h) PFER, q=1", "(i) PFER, q=5"
      )
    labeled <- mapply(function(g, lb) {
      ggdraw() +
        draw_grob(g) +
        draw_label(lb, x = 0.02, y = 0.98, hjust = 0, vjust = 1,
                   fontface = "bold", size = 22, color = "black")
    }, grobs, labels, SIMPLIFY = FALSE)
    
    combined <- plot_grid(plotlist = labeled, ncol = 3)
    ggsave(
      filename = file.path(save_dir, paste0(fig_title, "_combined.png")),
      plot = combined,
      width = 16, height = 16, dpi = 300, bg = "white"
    )
}