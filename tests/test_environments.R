## Test and learn scoping and environments
## I need to refine this and add more notes related to whay I learned

run.test.environments <- FALSE

if(run.test.environments){

  ## Scoping, frames and environments
  g <- function(a, b, c = 0) a + b + c
  print(parent.frame())  
  f1 <- function(){
    x <- 1
    print(ls())
    print(ls(envir = parent.frame()))
    print(parent.frame())
    f2 <- function(){
      w <- 2
      print(ls())
      print(ls(envir = parent.frame()))
      print(parent.frame())
      ans1 <- g(a = w, eval(x))
      f3 <- function(){
        print(ls(envir = parent.frame()))  
        print(parent.frame())
        print(eval(x))
        z <- eval(ans1) + eval(w) + eval(x)
      }
      ans1 <- f3()
    }
    ans <- f2()
    ans
  }
  
}
