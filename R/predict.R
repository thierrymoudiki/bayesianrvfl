predict_myridge <- function(obj, newx)
{
  my_scale(x = newx, xm = obj$xm, xsd = obj$scales)%*%obj$coef + obj$ym
}
