#
# Color keys used for plots and charts
#

# should share the same intensity to make them comparable
#
# Refs:
# http://colorbrewer2.org/
# http://www.google.com/url?sa=t&rct=j&q=chart%20color%20schemes&source=web&cd=2&cad=rja&ved=0CDcQxQEwAQ&url=http%3A%2F%2Fdocs.google.com%2Fviewer%3Fa%3Dv%26q%3Dcache%3Abx-BAPA77LcJ%3Awww.perceptualedge.com%2Farticles%2Fvisual_business_intelligence%2Frules_for_using_color.pdf%2Bchart%2Bcolor%2Bschemes%26hl%3Den%26gl%3Dus%26pid%3Dbl%26srcid%3DADGEESgPP_1I0BVOOhY3cGC3RxJ3nezMENv48fx2ltl-kgL2aQr6WfhaUZSVGGqzIz8KN5HgC8WPghpdK7fEybcKJ6-nYA-q9YNacjrjNBHKemSjD8s0zv64avlL88O5CKzPlXFuI9d_%26sig%3DAHIEtbSRp8IxbeRkwSrrhoEyGC_zoFkLtQ&ei=GBfvUPXADu6P0QHr7YGYDg&usg=AFQjCNF8y7mmVlL7ZKH1UEeHUFoOd4irdA

#' Common color scheme
#'
#' \describe{
#' \item{GROcap}{green}
#' \item{CAGE}{orange}
#' \item{GROseq}{purple}
#' \item{Fwd}{red}
#' \item{Rev}{blue}
#' }
#'
#' @export
gcol = list(
  grocap = "#4DAF4A",
  cage = "#FF7F00",
  groseq = "#984EA3",
  fwd = "#E41A1C",
  rev = "#377EB8")

#' Color ramp (between fwd and rev)
#'
#' @param n number of colors
#' @export
gcol.ramp <- function(n) {
  colorRampPalette(c(gcol$rev, gcol$fwd), space="rgb")(n)
}
