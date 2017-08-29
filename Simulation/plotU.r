MO <- 3
cf <- 2.5
x1 <- seq(0, 4, .05)
y1 <- -(x1-2)^2 + 4
y 	<- y1
sel <- x1 > MO
y[sel]	<- y1[sel] + cf
x	<- c(x1, MO)
y	<- c(y, -(MO-2)^2 +4 + cf)
ggtmp	<- data.frame(y = x, U = y)
ggplot(ggtmp, aes(y, U)) + geom_line()

alpha <- .2
plot(x1[sel], -alpha/log(x1[sel]-MO+1), type = "l")
y2	  <- y1
y2[sel] <- y2[sel] + cf - alpha/log(x1[sel] - MO + 1.1)
plot(x1, y2, type = "l")
x2		<- c(x1, MO)
y2		<- c(y2, -(MO-2)^2 +4 + cf - alpha/log(1.1))
ggtmp1	<- data.frame(y = x2, U = y2)
ggplot(ggtmp1, aes(y, U)) + geom_line() + 
	geom_line(data = subset(ggtmp, y>=MO), aes(y, U), linetype = 2) + 
	geom_vline(xintercept = MO, color = "red") + 
	scale_x_continuous(breaks = MO, labels = "MO") + 
	theme_bw() + 
	theme(axis.text.y = element_blank()) 
	