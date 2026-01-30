### script_intro_plot.R ---
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Jan  8 2026 (09:36) 
## Version: 
## Last-Updated: jan 28 2026 (08:47) 
##           By: Thomas Alexander Gerds
##     Update #: 196
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(data.table)
library(ggplot2)
library(ranger)
simulate_data <- function(n) {
    x <- runif(n)
    y <- rbinom(n = n, size = 1, prob = plogis(x^2))
    data <- data.table(x = x, y = factor(y,levels = c(0,1)))
    return(data)
}
set.seed(18)
d <- simulate_data(n = 1000)
set.seed(8)
## Fit correctly specified logistic regression model
fit_glm <- glm(y ~ I(x^2),data = d_boot,family = "binomial")
## Fit overfitting random forest model
fit_rf <- ranger::ranger(formula = y ~ x,
                         data = d_boot,
                         probability = TRUE,
                         num.trees = 1000,
                         mtry = 1,
                         min.node.size = 300)
## Non-parametric bootstrap with replacement
d_boot <- d[sample(1:nrow(d), size = nrow(d), replace = TRUE)]
## 10-fold CV on bootstrap sample
score_npboot <- riskRegression::Score(
                                    list("GLM" = fit_glm,
                                         "RF" = fit_rf),
                                    formula = y ~ 1,
                                    data = d_boot,
                                    contrasts = NULL,
                                    se.fit = FALSE,
                                    metrics = c("Brier", "AUC"),
                                    split.method = "cv10"
                                )
## Subsampling bootstrap 
d_subboot <- d[sample(1:nrow(d), size = 0.632*nrow(d), replace = FALSE)]
## 10-fold CV on bootstrap sample
score_subboot <- riskRegression::Score(
                                     list("GLM" = fit_glm,
                                          "RF" = fit_rf),
                                     formula = y ~ 1,
                                     data = d_subboot,
                                     contrasts = NULL,
                                     se.fit = FALSE,
                                     metrics = c("Brier", "AUC"),
                                     split.method = "cv10"
                                 )
summary(score_npboot)
summary(score_subboot)

set.seed(17)
x <- seq(0,1,0.001)
pred_df <- data.table(x = rep(x,2),
                      Model = rep(c("Data generating model","Random forest (n=1000)"),c(length(x),length(x))),
                      prediction = c(plogis(x^2),
                                     predict(fit_rf, data = data.table(x = x))$predictions[,2]))
g <- ggplot(pred_df, aes(x = x, y = prediction, color = Model,group = Model)) +
    geom_line(linewidth = 1) +
    labs(y =expression(paste("P(Y=1|X=x)=expit(",x^2,")")))+
    theme_bw() +
    theme(legend.position = c(.35,.85)) + theme(legend.title=element_blank())+
    theme(panel.grid = element_blank())+     theme(text = element_text(size=20))+
    scale_y_continuous(labels = scales::percent, limits = c(0,1))+
    scale_x_continuous(limits = c(0, 1)) +         
    scale_color_manual(values = c("#000000","#E69F00"),labels = c("Data generating model", "Random forest (n=1000)"))+
    theme(legend.key.height = unit(1.5, "cm"))
ggsave(g,filename = "./paper_statistics_in_medicine/figures/intro_plot.pdf", width = 9, height = 7, device = cairo_pdf)
g

### script_intro_plot.R ends here
