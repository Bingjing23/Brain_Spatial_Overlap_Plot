# Flatten brainMapR subcortical coordinates for the two standard views.
# The ROI translations follow brainMapR::addFlatCoordinatesSubcortical
# (Baptiste Couvy-Duchesne, MIT license).

add_flat_coordinates_subcortical <- function(annotation, hemi) {
  if (!hemi %in% c("lh", "rh")) stop("Unknown hemisphere: ", hemi)
  required <- c("ROINb", "X", "Y", "Z")
  missing <- setdiff(required, names(annotation))
  if (length(missing)) stop("Missing subcortical atlas columns: ", paste(missing, collapse = ", "))

  x <- data.table::copy(data.table::as.data.table(annotation))
  x[, `:=`(Xf = X, Yf = Y, Zf = Z)]
  hippocampus <- c(17L, 53L)
  amygdala <- c(18L, 54L)
  thalamus <- c(10L, 49L)
  caudate <- c(11L, 50L)
  accumbens <- c(26L, 58L)
  pallidum <- c(13L, 52L)

  if (hemi == "rh") {
    x[ROINb %in% hippocampus, Yf := Yf - 22]
    x[ROINb %in% amygdala, `:=`(Yf = Yf - 22, Zf = Zf + 9)]
    x[ROINb %in% thalamus, Zf := Zf - 22]
    x[ROINb %in% caudate, Yf := Yf + 24]
    x[ROINb %in% accumbens, Zf := Zf + 14]
    x[ROINb %in% pallidum, Yf := Yf - 17]
  } else {
    x[ROINb %in% hippocampus, Yf := Yf - 25]
    x[ROINb %in% amygdala, `:=`(Yf = Yf - 25, Zf = Zf + 9)]
    x[ROINb %in% thalamus, Zf := Zf - 26]
    x[ROINb %in% caudate, Yf := Yf + 14]
    x[ROINb %in% accumbens, Zf := Zf + 21]
    x[ROINb %in% pallidum, Yf := Yf - 24]
  }

  x[, `:=`(Xf2 = -Xf, Yf2 = Yf, Zf2 = -Zf)]
  if (hemi == "rh") {
    x[ROINb %in% accumbens, Zf2 := Zf2 - 7]
    x[ROINb %in% pallidum, Yf2 := Yf2 - 8]
    x[ROINb %in% amygdala, Yf2 := Yf2 - 5]
    x[ROINb %in% hippocampus, Yf2 := Yf2 - 5]
    x[ROINb %in% caudate, Yf2 := Yf2 - 8]
  } else {
    x[ROINb %in% caudate, Yf2 := Yf2 + 8]
    x[ROINb %in% pallidum, Yf2 := Yf2 + 6]
    x[ROINb %in% amygdala, Yf2 := Yf2 + 3]
    x[ROINb %in% hippocampus, Yf2 := Yf2 + 3]
    x[ROINb %in% accumbens, Zf2 := Zf2 + 6]
    x[ROINb %in% thalamus, Zf2 := Zf2 + 2]
  }
  x
}
