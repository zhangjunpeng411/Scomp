## Identifying ncRNA synergistic competing network based on ncRNA related ceRNA network by using sensitivity correlation
Scomp_sc <- function(ceRnet, ncRExpData, mRExpData, minSharedmR = 1, poscorcutoff = 0, pvaluecutoff = 0.05,
    senscorcutoff = 0.1) {

    m1 <- nrow(ceRnet)
    n1 <- ncol(ceRnet)    
    
    ncRExpDataNames <- as.matrix(colnames(ncRExpData))
    mRExpDataNames <- as.matrix(colnames(mRExpData))

    ncR <- ceRnet[, 1]
    mR <- ceRnet[, 2]

    ncRSym <- unique(ncR)
    mRSym <- unique(mR)

    m2 <- length(ncRSym)

    # Initialize variables
    ncRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
    C <- matrix(NA, m2 * (m2 - 1)/2, 5)

    for (i in seq_len(m2 - 1)) {
        for (j in seq(i + 1, m2)) {

            Interin1 <- ceRnet[which(ceRnet[, 1] %in% ncRSym[i]), 2]
            Interin2 <- ceRnet[which(ceRnet[, 1] %in% ncRSym[j]), 2]

            M1 <- length(Interin1)
            M2 <- length(Interin2)
            M3 <- length(intersect(Interin1, Interin2))
            M4 <- length(mRSym)
            M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)

            if (M3 >= minSharedmR & M5 < pvaluecutoff) {

                ncRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- ncRSym[i]
                ncRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- ncRSym[j]

                ncRExpIdx1 <- which(ncRExpDataNames %in% ncRSym[i])
                ncRExpIdx2 <- which(ncRExpDataNames %in% ncRSym[j])
		mRExpIdx <- which(mRExpDataNames %in% intersect(Interin1, Interin2))

                # Calculate sensitivity correlation of each ncRNA-ncRNA pair
                M6 <- cor.test(ncRExpData[, ncRExpIdx1], ncRExpData[, ncRExpIdx2])$estimate
                M7 <- cor.test(ncRExpData[, ncRExpIdx1], ncRExpData[, ncRExpIdx2])$p.value
		M8 <- M6 - pcor.shrink(cbind(ncRExpData[, ncRExpIdx1], ncRExpData[, ncRExpIdx2],
                  mRExpData[, mRExpIdx]), verbose = FALSE)[1, 2]

                C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M3
                C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M5
                C[(i - 1) * m2 + j - sum(seq_len(i)), 3] <- M6
                C[(i - 1) * m2 + j - sum(seq_len(i)), 4] <- M7
                C[(i - 1) * m2 + j - sum(seq_len(i)), 5] <- M8
            }
        }
    }

    # Extract ncRNA-ncRNA pairs with sensitivity correlation more than senscorcutoff.
    ncRInt <- ncRInt[which((C[, 2] < pvaluecutoff & C[, 3] > 
        poscorcutoff & C[, 4] < pvaluecutoff & C[, 5] > senscorcutoff) == "TRUE"), ]

    C <- C[which((C[, 2] < pvaluecutoff & C[, 3] > 
        poscorcutoff & C[, 4] < pvaluecutoff & C[, 5] > senscorcutoff) == "TRUE"), ]

    if (is.vector(C)) {
        SPPCncRInt <- c(ncRInt, C)
        names(SPPCceRInt) <- c("ncRNA_1", "ncRNA_2", "#shared mRNAs", "pvalue of shared mRNAs",
            "correlation", "pvalue of positive correlation", "sensitivity pearson correlation")
    } else {
        SPPCncRInt <- cbind(ncRInt, C)
        colnames(SPPCncRInt) <- c("ncRNA_1", "ncRNA_2", "#shared mRNAs", "pvalue of shared mRNAs",
            "correlation", "pvalue of positive correlation", "sensitivity pearson correlation")
    }

    return(SPPCncRInt)
}

Random_net_clusterCoeff <- function(nodes.num, edges.num, perm = 10000, directed = FALSE) {
    set.seed(123)
    res <- c()
    for (i in seq(perm)) {
    g <- sample_pa(n = nodes.num, m = edges.num, directed = directed)
    g <- delete_edges(g, sample(1:gsize(g), size = gsize(g) - edges.num))
        res[i] <- transitivity(g, type="average")
    }

    return(list(mean(res), sd(res)))
}

## Survival analysis of modules
moduleSurvival <- function(Modulelist, ExpData, SurvData, devidePercentage = 0.5, plot = FALSE) {

    ExpDataNames <- colnames(ExpData)    
    myfit <- list()
    LogRank <- list()

    for (i in seq_along(Modulelist)) {
        Interin_Data <- cbind(SurvData[, seq(2, 3)], ExpData[, which(ExpDataNames %in% Modulelist[[i]])])
        Interin_Data <- na.omit(Interin_Data)

        try_mm <- try(coxph(survival::Surv(time, status) ~ ., data = data.frame(Interin_Data)),
            silent = TRUE)
        if ("try-error" %in% class(try_mm))
            next

        mm <- coxph(survival::Surv(time, status) ~ ., data = data.frame(Interin_Data))

        Risk_score <- predict(mm, newdata = data.frame(Interin_Data), type = "risk")

        group <- rep("NA", dim(Interin_Data)[1])
        group[Risk_score > quantile(Risk_score, probs = devidePercentage)] <- "High"
        group[Risk_score <= quantile(Risk_score, probs = devidePercentage)] <- "Low"

        Data <- cbind(Interin_Data[, seq_len(2)], group)
        myfit[[i]] <- survfit(survival::Surv(time, status) ~ group, data = Data)

        sdf <- survdiff(survival::Surv(time, status) ~ group, data = Data)
        sdf.p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
        HR <- (sdf$obs[1]/sdf$exp[1])/(sdf$obs[2]/sdf$exp[2])
        HRlow95 <- exp(log(HR) - qnorm(0.975) * sqrt(1/sdf$exp[1] + 1/sdf$exp[2]))
        HRup95 <- exp(log(HR) + qnorm(0.975) * sqrt(1/sdf$exp[1] + 1/sdf$exp[2]))

        LogRank[[i]] <- c(sdf$chisq, sdf.p.val, HR, HRlow95, HRup95)
    }

    if (plot) {
        for (i in seq_along(myfit)) {
            if (!is.null(LogRank[[i]])) {
                dev.new()
                plot(myfit[[i]], lty = 1, col = c("red", "green"), main = paste("Module", i), xlab = "Time (Months)",
                    ylab = "Probability of survival")

                legend("bottomleft", legend = c("High risk group", "Low risk group"), lty = seq_len(2),
                    col = c("red", "green"), bty = "n")
            }
        }
    }

    LogRank_res <- do.call(rbind, LogRank)

    if (length(myfit) >= 1) {
        colnames(LogRank_res) <- c("Chi-square", "p-value", "HR", "HRlow95", "HRup95")
        names(LogRank) <- seq_along(Modulelist)
        LogRank[sapply(LogRank, is.null)] <- NULL
        rownames(LogRank_res) <- paste("Module", names(LogRank))
    }

    return(LogRank_res)
}


## Evaluate the performance of each module for classifying tumor types
module.classify <- function(Exp, tumor_type, Modulelist, method = "br", base.algorith = "SVM", cv.folds = 10, 
	                    cv.sampling = "stratified", cv.seed = 12345) {

    module_Exp <- lapply(seq_along(Modulelist), function(i) Exp[, which(colnames(Exp) %in% Modulelist[[i]])])
    unique_type <- unique(tumor_type[, 2])
    class_infor <- do.call(cbind, lapply(seq(unique_type), function(i) as.numeric(tumor_type[, 2] == unique_type[i])))
    module_classify <- list()

    for (i in seq_along(Modulelist)){        
	
        temp <- as.data.frame(cbind(module_Exp[[i]], class_infor))
	Indices <- ncol(temp)
	temp_mldr <- mldr_from_dataframe(temp, labelIndices = (Indices-length(unique_type)+1):Indices, name = "TEMPMLDR")
        temp_res <- cv(temp_mldr, method = method, base.algorith = base.algorith, cv.folds = cv.folds, 
	               cv.sampling = cv.sampling, cv.seed = cv.seed)
        module_classify[[i]] <- temp_res

    }

    return(module_classify)
}

