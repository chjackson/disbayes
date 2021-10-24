test_that("ci2num", {
    
    est <- c(0.2, 0.3, 0.4)
    lower <- c(0.1, 0.2, 0.3)
    upper <- c(0.3, 0.4, 0.5)
    expect_equal(ci2num(est, lower, upper)$num, c(13,24,36))
    
    ciz1 <- ci2num(0, 0, 1, denom0=100)
    expect_equal(ciz1$num, 0)
    expect_equal(ciz1$denom, 100)
    
    expect_equal(ci2num(0.1, 0, 1)$num, 9)
    
    expect_error(ci2num(-1, 0.1, 0.2), "should be in \\[0,1\\]")
    expect_error(ci2num(0.1, 2, 0.2), "should be in \\[0,1\\]")
    expect_error(ci2num(0.1, 0.01, 2), "should be in \\[0,1\\]")
    expect_error(ci2num(0.01, 0.1, 0.2), "should be inside the credible interval")
    expect_error(ci2num(0.01, 0.2, 0.1), "should be < upper")
})    
