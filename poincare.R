
library("pracma")
S_T <- function(SH, SD) {
    cross(SH, SD)
}

library("rgl")
require("sphereplot")

set.seed(101)

x = 1.5
y = 0
z = 0

#' @export
rgl.sphere.grid <- function(radius = 1, col.long = "blue", col.lat = "blue",
                           deggap = 15, textgap = 90, longtype = "D",
                           add = FALSE, radaxis = FALSE,
                           radlab = "Radius", gridtype = "wire") {
    ## This code is adapted from the rgl.sphgrid function in
    ## package sphereplot by Aaron Robotham

    for (lat in seq(-90, 90, by = deggap)) {
        if (lat == 0) {
            col.grid = "grey50"
        }
        else {
            col.grid = "grey"
        }
        rgl::plot3d(sphereplot::sph2car(long = seq(0, 360, len = 100), lat = lat,
                                    radius = radius, deg = T), col = col.grid, add = T, type = "l")
    }
    for (long in seq(0, 360 - deggap, by = deggap)) {
        if (long == 0) {
            col.grid = "grey50"
        }
        else {
            col.grid = "grey"
        }
        rgl::plot3d(sphereplot::sph2car(long = long, lat = seq(-90, 90, len = 100),
                                    radius = radius, deg = T), col = col.grid, add = T, type = "l")
    }
    if (longtype == "H") {
        scale = 15
    }
    if (longtype == "D") {
        scale = 1
    }
    sphereplot::rgl.sphtext(long = 0, lat = seq(-90, 90, by = textgap * 2),
                          radius = 1.1 * radius, text = c("|L>", "|R>"), #seq(-90, 90, by = textgap*2),
                          deg = TRUE, col = col.lat)

    sphereplot::rgl.sphtext(long = seq(0, 360 - textgap, by = textgap), lat = 0,
                          radius = 1.1 * radius, text = c("|H>", "|D>", "|V>", "|A>"), #seq(0, 360 - textgap,
                                                          #by = textgap)/scale,
  deg = TRUE, col = col.long)

    if (radaxis) {

        radpretty = pretty(c(0, radius))
        radpretty = radpretty[radpretty <= radius]

        rgl::lines3d(c(0, 0), c(0, max(radpretty)), c(0, 0), col = "grey50")

        for (i in 1:length(radpretty)) {
            rgl::lines3d(c(0, 0), c(radpretty[i], radpretty[i]),
                   c(0, 0, radius / 50), col = "grey50")
            rgl::text3d(0, radpretty[i], radius / 15, radpretty[i],
                  col = "darkgreen")
        }

        rgl::text3d(0, radius / 2, - radius / 25, radlab)
    }
}

retarder <- function(retard_, phi_) {
    C2 <- cos(2 * phi_)#+pi)
    S2 <- sin(2 * phi_)#+pi)  # + pi is added due to angle reference to be horizontal
    M <- matrix(nrow = 3, ncol = 3)

    M[1, 1] <- C2 * C2 + S2 * S2 * cos(retard_)
    M[1, 2] <- S2 * C2 * (1 - cos(retard_))
    M[1, 3] <- S2 * sin(retard_)


    M[2, 1] <- S2 * C2 * (1 - cos(retard_))
    M[2, 2] <- S2 * S2 + C2 * C2 * cos(retard_)
    M[2, 3] <- -C2 * sin(retard_)


    M[3, 1] <- -S2 * sin(retard_)
    M[3, 2] <- C2 * sin(retard_)
    M[3, 3] <- cos(retard_)

    M

}

cartesian_to_spherical <- function(S1, S2, S3) {
    # return degree angle

    r = sqrt(S1 * S1 + S2 * S2 + S3 * S3)
    chi2 = asin(S3 / r) * 180 / pi
    #if(chi2 != 0)
        #psi2 = asin(S2 / (r * cos(chi2 * pi / 180))) * 180 / pi
    #else
        #psi2 = 0
    psi2 = atan(S2 / S1) * 180 / pi

    if (S1 <= 0 && S2 < 0){
      return(c(r, chi2, psi2 - 180))
    }  # 3rd quadrant

    else if (S1 < 0 && S2 >= 0){
      return(c(r, chi2, psi2 + 180))
    }  # 2nd quadrant

    else if (S1 == 0 && S2 > 0){
      return(c(r, chi2, 90))
    }

    else if (S1 == 0 && S2 < 0){
      return(c(r,chi2,-90))
    }

    else{
      return(c(r,chi2,psi2))
    }


}

spherical_to_cartesian <- function(r, ellipticity, azimuth) {
    #S=c(S1,S2,S3)
    S3 = r * sin(ellipticity * pi / 180)
    S1 = r * cos(ellipticity * pi / 180) * cos(azimuth * pi / 180)
    S2 = r * cos(ellipticity * pi / 180) * sin(azimuth * pi / 180)
    return(c(S1, S2, S3))
}

SOP <- function ( sph_SOP, color){
  plot3d(sph_SOP[1]*cos(sph_SOP[2] * pi / 180) * cos(sph_SOP[3] * pi / 180), sph_SOP[1]*cos(sph_SOP[2] * pi / 180) * sin(sph_SOP[3] * pi / 180), sph_SOP[1]*sin(sph_SOP[2] * pi / 180), xlab = '', ylab = '', zlab = '', type = 'p', size = 10, col = color, lwd = 2, xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), axes = F)
}

color.gradient <- function(x, colors = c("red", "yellow"), colsteps = 100) {
    return(colorRampPalette(colors)(colsteps)[findInterval(x, seq(min(x), max(x), length.out = colsteps))])
}


find_Stokes_vector <- function(DOP, ellipticity, azimuth){
  S_vector = spherical_to_cartesian(DOP, ellipticity, azimuth)
  return(S_vector)
}

poincare <- function(chi2_init, psi2_init, retard) {
    # input angle is in degree

    df <- data.frame(S1 = double(), S2 = double(), S3 = double())
    colfunc <- colorRampPalette(c("red", "blue"))
    # initial Stokes vector
    S_init <- c( cos(chi2_init * pi / 180) * cos(psi2_init * pi / 180), sin(psi2_init * pi / 180) * cos(chi2_init * pi / 180), sin(chi2_init * pi / 180))

    #delta_in_degree = 30*i # retardation in degree
    for (i in 0:179) {
        phi <- i

        # final Stokes vector
        S_final <- retarder(retard * pi / 180, phi * pi / 180) %*% S_init

        df <- rbind(df, data.frame(S_final[1, 1], S_final[2, 1], S_final[3, 1]))
    }
    {
        library(scatterplot3d)
        #    par(bg = NA)
        library(rgl)
        require(rgl)
        rgl.material(front = "points", back = "points")
        plot3d(df, type = 'l', col = color.gradient(as.numeric(rownames(df))), lwd = 5, size = 5, xlab = '', ylab = '', zlab = '', xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), axes = F)
    }
}



SH <- c(1,2.077*2,-0.579*2)
SD <- c(1,0.371*2,45.9*2)
SV <- c(1,0.172*2,87.984*2 )
SA <- c(1,-1.239*2,-43.85*2 )
SOP(SD,'red')
SOP(SH, 'red')
SOP(SV, 'red')
SOP(SA, 'red')

rgl.sphere.grid()
