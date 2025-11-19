######################## source common libraries ###############################
source("common.R")

######################## import input ##########################################
# this study's patients
data <- fread("../result/Fig2/num_patients_per_biological_category.csv")

# this study's patients, classified as per Euroclass
euroclass_STUDYID_df <- fread("../result/euroclass_STUDYID.csv")

# Euroclass patients count
euroclass_patient_count_df <- fread("../input/Wehr_al_2028_patient_count.csv")

############# prep ############################################################

data$study <- "INTREPID"

data <- data %>% 
  dplyr::select(biological_level, number_of_patients, study, biological_measure)

# include Euroclass categories and Euroclass counts
edf <- euroclass_STUDYID_df %>%
  mutate(
    euroclass_group = case_match(
      euroclass_group,
      "smBminusTrhi_Euroclass"    ~ "smB-Trhi",
      "smBminusTrnorm_Euroclass"  ~ "smB-Trnorm",
      "smBplus21lo_Euroclass"     ~ "smB+CD21lo",
      "smBplus21norm_Euroclass"   ~ "smB+CD21norm",
      "smBminusCD21lo_Euroclass"  ~ "smB-CD21lo",
     "smBminusCD21norm_Euroclass" ~ "smB-CD21norm",
   .default = euroclass_group)) %>%
  dplyr::filter(euroclass_group %in% c("smB-Trhi",
                                       "smB-Trnorm",
                                       "smB+CD21lo",
                                       "smB+CD21norm",
                                       "smB-CD21lo",
                                       "smB-CD21norm"))



for(bl in euroclass_patient_count_df$biological_level){
  n <- euroclass_patient_count_df %>% 
    dplyr::filter(biological_level == bl) %>%
    dplyr::select(number_of_patients) %>% unlist()
  data <- data %>% 
    add_row(biological_level = bl,
            number_of_patients = n,
            study ="EUROClass",
            biological_measure = bl)
}

euroclass_markers <- unique(edf$euroclass_group)

for(bl in euroclass_markers){
  n <- edf %>% 
    dplyr::filter(euroclass_group == bl) %>%
    nrow()
  data <- data %>% 
    add_row(biological_level = bl, 
            number_of_patients = n,
            study = "INTREPID",
            biological_measure = bl)
}


############# Fig2 ############################################################
df2 <- data %>% 
  dplyr::filter(biological_measure %in% c("CD21", "SmB", "Tr")) %>% 
  group_by(biological_measure) %>% 
  mutate(groupPerc = (number_of_patients / sum(number_of_patients)) * 100) %>%
  ungroup()

Fig2A <- df2 %>% 
  ggplot(aes(fill = biological_level,
             y = number_of_patients,
             x = biological_measure))  + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(
    label = paste0(
      round(groupPerc, 1),
      "%")),
    position = position_dodge(0.9),
    size = 3,
    vjust = -0.5 )+
  ylab("Number of patients") +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.text = element_text(size=8)) +
  scale_fill_manual(name = "",
                    values = Bcellcol,
                    labels = c("expansion of low CD21" = expression(CD21^norm),
                               "CD21 normal" = expression(CD21^low),
                               "SmB minus" = "smB-",
                               "SmB plus" = "smB+",
                               "Transitional B normal" = expression(Tr^norm),
                               "Transitional B high" = expression(Tr^hi))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  guides(fill = guide_legend(nrow = 2))




############# Fig2B ############################################################


# 3 B cell groups
Bcell_groups <- list(
  "smB+CD21" = c("smB+CD21lo", "smB+CD21norm"),
  "smB-CD21" = c("smB-CD21lo", "smB-CD21norm"),
  "smB-Tr"   = c("smB-Trnorm", "smB-Trhi")
)


Fig2Bplot_data <- data %>%
  dplyr::filter(biological_measure %in% unlist(Bcell_groups)) %>%
  mutate(Bcell_group = case_when(
      biological_measure %in% Bcell_groups[["smB+CD21"]] ~ "smB+CD21",
      biological_measure %in% Bcell_groups[["smB-CD21"]] ~ "smB-CD21",
      biological_measure %in% Bcell_groups[["smB-Tr"]]   ~ "smB-Tr")) %>%
  group_by(Bcell_group, study) %>%
  mutate(percentage_within_study =
           number_of_patients / sum(number_of_patients) * 100) %>%
  ungroup() %>%
  mutate(x_axis_label = paste(Bcell_group, study, sep = "\n")) %>%
  mutate(x_axis_label = factor(x_axis_label,
                               levels = c(
      "smB+CD21\nINTREPID", "smB+CD21\nEUROClass",
      "smB-CD21\nINTREPID", "smB-CD21\nEUROClass",
      "smB-Tr\nINTREPID", "smB-Tr\nEUROClass")))

smBTrpv <- Fig2Bplot_data %>%
  dplyr::filter(Bcell_group == "smB-Tr") %>%
  summarise(
    p_value = {ct <- xtabs(number_of_patients ~ study + biological_level,
                           data = .)
    chisq.test(ct)$p.value}) %>% as.numeric()

smBCD21pv <-Fig2Bplot_data %>%
  dplyr::filter(Bcell_group == "smB-CD21") %>%
  summarise(
    p_value = {ct <- xtabs(number_of_patients ~ study + biological_level,
                           data = .)
    chisq.test(ct)$p.value})%>% as.numeric()

smBpCD21pv <-Fig2Bplot_data %>%
  dplyr::filter(Bcell_group == "smB+CD21") %>%
  summarise(
    p_value = {ct <- xtabs(number_of_patients ~ study + biological_level,
                           data = .)
    chisq.test(ct)$p.value})%>% as.numeric()


stats_summary <- tibble(Bcell_group=c("smB+CD21", "smB-CD21","smB-Tr"),  
                        p_value= c(smBpCD21pv,smBCD21pv,smBTrpv),
                        y.position=c(95, 91, 146),
                        xmin = c("smB+CD21\nINTREPID",
                                 "smB-CD21\nINTREPID",
                                 "smB-Tr\nINTREPID"),
                        xmax = c("smB+CD21\nEUROClass",
                                 "smB-CD21\nEUROClass",
                                 "smB-Tr\nEUROClass"),
                        label = paste("\u03C7\u00B2 test, p =",
                                      round(p_value, digits = 3)))
  
  


all_colors <- c("smB+CD21lo"   = smBPLUSCD21col[1],
  "smB+CD21norm" = smBPLUSCD21col[2],
  "smB-CD21lo"   = smBMINUSCD21col[1],
  "smB-CD21norm" = smBMINUSCD21col[2],
  "smB-Trnorm"   = smBMINUSTrcol[1],
  "smB-Trhi"     = smBMINUSTrcol[2])

all_labels <- c("smB+CD21lo"   = expression("smB+CD21"^lo),
              "smB+CD21norm" = expression("smB+CD21"^norm),
              "smB-CD21lo"   = expression("smB-CD21"^lo),
              "smB-CD21norm" = expression("smB-CD21"^norm),
              "smB-Trnorm"   = expression("smB-Tr"^norm),
              "smB-Trhi"     = expression("smB-Tr"^hi))


Fig2B <- ggplot(Fig2Bplot_data,
       aes(x = x_axis_label,
           y = number_of_patients,
           fill = biological_level)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", percentage_within_study)),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3) +
  geom_bracket(data = stats_summary,
    inherit.aes = FALSE,
    aes(xmin = xmin,
        xmax = xmax,
        y.position = y.position,
        label = label),
    label.size = 3) +
  scale_fill_manual(values = all_colors,
                    labels = all_labels,
                    name = "",
                    breaks = names(all_colors)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(labels =
                     c("smB+CD21\nINTREPID" = "INTREPID",
                      "smB-CD21\nINTREPID"  = "INTREPID",
                      "smB-Tr\nINTREPID"    = "INTREPID",
                      "smB+CD21\nEUROClass" = "EUROClass",
                      "smB-CD21\nEUROClass" = "EUROClass",
                      "smB-Tr\nEUROClass"   = "EUROClass")) +
  labs(y = "Number of patients") +
  theme(axis.title.x   = element_blank(),
        legend.position = "bottom",
        legend.text     = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 2)) 


############# Layout ###########################################################

toppatch <- Fig2A
bottompatch <- Fig2B


patchwork1 <-  toppatch / bottompatch

patchwork1 + plot_annotation(tag_levels = 'a')

ggsave("../result/Fig2/Fig2.jpeg",
       width = 15,
       height = 15,
       units = "cm")

############# Legend details ###################################################
#c hsqr
stats_summary$Bcell_group[stats_summary$p_value<0.05]
