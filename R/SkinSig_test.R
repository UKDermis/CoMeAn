

SkinSig_annotation <- read.csv("SkinSig_annotation.csv")
unq <- unique(SkinSig_annotation$SkinSig.signature)
colors <- c("black", "#68b2f7", "orange", "#eda65f",
                  "#46f2a2", "#c76e20", "#d12a6d", "#d41f15",
                  "darkblue", "#79ab22", "pink", "#e647d5",
                  "#821522", "magenta", "#f21679", "#81eb3b",
                  "#f5a6d4", "#d1b0c3", "#f0210e", "white")
print(unq)
print(length(unq))
print(length(colors))

names(colors) <- unq

print(colors[[unq[1]]])
print(SkinSig_annotation$SkinSig.signature[1])
print(colors[[SkinSig_annotation$SkinSig.signature[1]]])

N <- nrow(SkinSig_annotation)
scols <- vector("list", N)
for(i in 1:N){
  scols[i] <- colors[[SkinSig_annotation$SkinSig.signature[i]]]
}
print(scols)
