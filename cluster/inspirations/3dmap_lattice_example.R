cloud(density ~ long + lat, state.info,
      subset = !(name %in% c("Alaska", "Hawaii")),
      panel.3d.cloud = function(...) {
          panel.3dmap(...)
          panel.3dscatter(...)
      }, 
      type = "h", scales = list(draw = FALSE), zoom = 1.1,
      xlim = state.map$range[1:2], ylim = state.map$range[3:4],
      xlab = NULL, ylab = NULL, zlab = NULL,
      aspect = c(diff(state.map$range[3:4]) / diff(state.map$range[1:2]), 0.3),
      panel.aspect = 0.75, lwd = 2, screen = list(z = 30, x = -60),
      par.settings = list(axis.line = list(col = "transparent"),
                          box.3d = list(col = "transparent", alpha = 0)))

