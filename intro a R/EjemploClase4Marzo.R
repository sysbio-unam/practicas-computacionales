> x <- 5
> y <- 16
> x+y

> x-y

> x*y

> y/x

> y%/%x

> y%%x

> y^x


> x<y

> x>y

> x<=5

> y>=20

> y == 16

> x != 5


> x <- c(2,8,3)
> y <- c(6,4,1)
> x+y

> x>y



> x <- c(TRUE,FALSE,0,6)
> y <- c(FALSE,TRUE,FALSE,TRUE)
> !x

> x&y

> x&&y

> x|y

> x||y


> x <- 5
> x

> x = 9
> x

> 10 -> x
> x


x <- -5
if(x > 0){
  print("Non-negative number")
} else {
  print("Negative number")
}


x <- runif(100,min=1,100)
count <- 0
for (val in x) {
  if(val == 0)  count = count+1
}
print(count)


i <- 1
while (i < 6) {
  print(i)
  i = i+1
}
