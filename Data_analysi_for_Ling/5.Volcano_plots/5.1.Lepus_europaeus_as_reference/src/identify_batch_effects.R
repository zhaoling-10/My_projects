# Make syntactically valid names for cross levels
metadata$cross <- make.names(metadata$cross)

# Model matrices for linear modeling
mod <- model.matrix(~0 + cross, data = metadata)  # Full model
mod0 <- model.matrix(~1, data = metadata)         # Null model

# Number of surrogate variables
n.sv <- num.sv(filtered_data, mod, method = "leek")
n.sv <-5
# Surrogate Variable Analysis (SVA)
svobj <- sva(as.matrix(filtered_data), mod, mod0, n.sv = n.sv)

# Model with surrogate variables
modSv <- cbind(mod, svobj$sv)
mod0Sv <- cbind(mod0, svobj$sv)

# Check the column names in the mod matrix
print(colnames(mod))

# Define contrast matrix using the actual column names from 'mod'
contrast.matrix <- limma::makeContrasts(
  tt8_4x_vs_2x_2x = `crossX2x_tt8.x.4x_Col.0` - `crossX2x_Col.0.x.2x_Col.0`,
  `2x_4x_vs_2x_2x` = `crossX2x_Col.0.x.4x_Col.0` - `crossX2x_Col.0.x.2x_Col.0`,
  `tt8_2x_vs_2x_2x` = `crossX2x_tt8.x.2x_Col.0` - `crossX2x_Col.0.x.2x_Col.0`,
  `4x_2x_vs_2x_2x` = `crossX4x_Col.0.x.2x_Col.0` - `crossX2x_Col.0.x.2x_Col.0`,
  `tt8_4x_vs_2x_4x` = `crossX2x_tt8.x.4x_Col.0` - `crossX2x_Col.0.x.4x_Col.0`,
  levels = mod  # Use 'mod' matrix without surrogate variables
)

# Linear model fitting using limma
fit <- limma::lmFit(as.matrix(filtered_data), modSv)

# Fit the contrast matrix and compute eBayes
fitContrasts <- limma::contrasts.fit(fit, contrast.matrix)
eb <- eBayes(fitContrasts)

# Get top table results
topTableF(eb, adjust = "BH")