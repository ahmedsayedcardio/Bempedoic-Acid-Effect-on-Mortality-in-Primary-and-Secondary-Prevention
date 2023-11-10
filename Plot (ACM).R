##Plot posterior distributions for all-cause mortality##

###Make draws of BA's treatment effect from the posterior distribution 
###under primary and secondary prevention settings
#Primary
primary_draws <- posterior_linpred(
  b_model_acm, #This is the model we used
  newdata = data.frame(ba = c(1, 0), #Calculate probabilities according to the use (or not) of BA
                       prevention_status = "Primary"), #In primary prevention
  transform = TRUE #Create predictions on the probability (rather than log-odds) scale
)
#Secondary (as above but change prevention status to "secondary")
secondary_draws <- posterior_linpred(
  b_model_acm, newdata = data.frame(ba = c(1, 0), prevention_status = "Secondary"),
  transform = TRUE
)
#The output is the probability (risk) of the outcome (death) with BA (1st column) and without BA (2nd column)


###Calculate the ratio of these two probabilities/risks (the risk difference; rd)
primary_draws <- primary_draws %>% data.frame %>% mutate(rd = .[, 1] - .[, 2], #Dividing one probability by the other yields the rd
                                                         prevention_status = "Primary") #Under primary prevention
secondary_draws <- secondary_draws %>% data.frame %>% mutate(rd = .[, 1] - .[, 2], #As above, this yields the rd
                                                             prevention_status = "Secondary") #Under secondary prevention
###We now have a distribution of risk differences under the primary and secondary prevention settings
###This is needed to for the plot of the distribution of BA's effect that we will create a little further below


###Now, we need to calculate the point estimate for the risk difference and the 95% upper and lower bounds
#For primary prevention
primary_rd <- posterior_summary(primary_draws %>% select(rd))
#For secondary prevention
secondary_rd <- posterior_summary(secondary_draws %>% select(rd))
#The output is the risk difference under primary and secondary prevention settings

#Finally, calculate the probability of benefit (reducing ACM) under the 1ry prevention setting
#and the probabiltiy of harm (increasing ACM) under the 2ry prevention setting
primary_post_prob <- b_model_acm %>% hypothesis("ba < 0") %>%
  {.$hypothesis$Post.Prob} %>% #Extract the posterior probability 
  {. * 100} %>% #Convert to % for ease of interpretation
  round(., 2) %>% format(., nsmall = 2) #Round to 2 decimal places (and show trailing zeros)
secondary_post_prob <- b_model_acm %>% hypothesis("ba + ba:prevention_statusSecondary > 0") %>%
  {.$hypothesis$Post.Prob} %>% #Extract the posterior probability 
  {. * 100} %>% #Convert to % for ease of interpretation
  round(., 2) %>% format(., nsmall = 2) #Round to 2 decimal places (and show trailing zeros)


#Calculate the probability that bempedoic acid has a less favorable effect in the secondary prevention setting
#(that is, the probability that it either increases mortality or has a lesser mortality benefit)
interaction_post_prob <- b_model_acm %>% hypothesis(hypothesis = "ba:prevention_statusSecondary > 0") %>%
  {.$hypothesis$Post.Prob} %>% #Extract the posterior probability 
  {. * 100} %>% #Convert to % for ease of interpretation
  round(., 2) %>% format(., nsmall = 2) #Round to 2 decimal places (and show trailing zeros)

##Plot the resulting posterior distributions for BA's effect on 1ry and 2ry prevention
#First, let's the primary and secondary prevention settings into a single dataframe
combined_draws <- bind_rows(primary_draws, secondary_draws)

#Give the plot a title and a subtitle (the subtitle is composed of two parts)
main_title <- "Effect of Bempedoic Acid on All-Cause Mortality in the Primary and Secondary Prevention Settings."
subtitle_p1 <- "Effects were obtained from a logistic regression model with death as a binary outcome and use of bempedoic acid, primary versus secondary prevention status, and an interaction between them as predictors."
subtitle_p2 <- paste0("\nThe posterior probability of a decrease in mortality in the primary prevention setting was ", primary_post_prob, "%. The posterior probability of an increase in mortality in the secondary prevention setting was ", secondary_post_prob, "%.\nThe posterior probability of a less favorable effect in the secondary prevention versus primary prevention setting (lesser benefit or harm) was ", interaction_post_prob, "%.")
subtitle <- paste0(subtitle_p1, subtitle_p2)

#Plot
acm_plot <- ggplot(data = combined_draws,
       aes(y = 1,
           x = rd*100,
           fill = prevention_status
       )) +
  #Add Density plots
  stat_halfeye(alpha = 0.75, .width = 0.95) +
  #Set colors
  scale_fill_tableau(name = "Prevention status") +
  #Create title
  ggtitle(main_title,
          subtitle = subtitle) +
  geom_vline(xintercept = 0, color = "black", 
             lwd = 0.75, linetype = 1) +
  #X and Y axes aesthetics
  scale_y_discrete(name = NULL, expand = c(0, 0.03)) +
  scale_x_continuous(name = "Risk Difference (%)",
                     expand = c(0, 0),
                     breaks = seq(-10, 10, 1)) +
  #Set theme
  theme_pubclean() +
  theme(text = element_text(size = 23),
        plot.title=element_text(face = "bold", hjust = 0.0, size = 18),
        plot.subtitle = element_text(face = "bold", size = 10, hjust = 0.0, color = "grey45"),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black", size = 1.2),
        plot.margin = margin(0.5, 1, 0.5, 1, "cm"),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "bottom",
        legend.text = element_text(size = 12, face = "bold"),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.75, "cm")
        )


ggsave(plot = acm_plot,
       filename = "C:/Ahmed's Stuff/ResearchStuff/BA_Primary_v_Secondary/Figure 1.png",
       dpi = 600,
       width = 16,
       height = 9)
