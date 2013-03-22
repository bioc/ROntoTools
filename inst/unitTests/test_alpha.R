test_alpha1MR <- function()
{
  
  checkEquals(alpha1MR(.5), 0)
  
  c <- (1:10)/10
  
  checkEquals(alpha1MR(c), 1-c/max(c))
  
  checkException(alpha1MR("c"), "pv is not numeric")
  
}

test_alphaMLG <- function()
{
  
  checkEquals(alphaMLG(.5), 0)
  
  c <- (1:10)/10
  
  checkEquals(alphaMLG(c), -log10(c/max(c)))
  
  checkException(alphaMLG("c"), "pv is not numeric")
  
}