aggregate_transmitted_f <- rowSums(false_pos_pf)
percentages_transmitted_f <- aggregate_transmitted_f/sum(aggregate_transmitted_f)
per_sample_false_pos_rate_f <- apply(false_pos_pf, 2, function(xx) xx[1]/sum(xx))

aggregate_transmitted_m <- rowSums(false_pos_pm)
percentages_transmitted_m <- aggregate_transmitted_m/sum(aggregate_transmitted_m)
per_sample_false_pos_rate_m <- apply(false_pos_pm, 2, function(xx) xx[1]/sum(xx))
