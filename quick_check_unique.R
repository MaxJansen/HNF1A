#This is a script to check duplicates in names within and between batches.
#Use in "2.large_table.R" script.
#I knew they weren't unique within batches, check if between batches:
Batch1_select <- Mastertable[(Mastertable$BatchName == "Batch1"), ]
Batch1_select <- Batch1_select[!duplicated(Batch1_select$Name), ]

Batch2_select <- Mastertable[(Mastertable$BatchName == "Batch2"), ]
Batch2_select <- Batch2_select[!duplicated(Batch2_select$Name), ]

Batch2second_select <- Mastertable[(Mastertable$BatchName == "Batch2second"), ]
Batch2second_select <- Batch2second_select[!duplicated(Batch2second_select$Name), ]
Batch2_new_select <- rbind(Batch2_select,Batch2second_select)
Batch2_new_select <- Batch2_new_select[!duplicated(Batch2_new_select$Name), ]

Batch3_select <- Mastertable[(Mastertable$BatchName == "Batch3"), ]
Batch3_select <- Batch3_select[!duplicated(Batch3_select$Name), ]
Batch4_select <- Mastertable[(Mastertable$BatchName == "Batch4"), ]
Batch4_select <- Batch4_select[!duplicated(Batch4_select$Name), ]

full_bind <- rbind(Batch1_select, Batch2_new_select, Batch3_select, Batch4_select)
ultimate_test <- full_bind[!duplicated(full_bind$Name), ]
Batch4_select <- Mastertable[(Mastertable$BatchName == "Batch4"), ]