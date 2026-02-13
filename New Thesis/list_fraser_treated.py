import pandas as pd
df = pd.read_csv("merged_full_data.csv")
df = df[(df.Year >= 1995) & (df.Year <= 2022) & df.ISO3.notna()]
nobs = df.groupby("ISO3").size()
panel = nobs[nobs >= 15].index.tolist()
df = df[df.ISO3.isin(panel)]
df = df.sort_values(["ISO3", "Year"])
df["lag_p"] = df.groupby("ISO3")["fraser_Summary"].shift(1)
df["jump"] = df["fraser_Summary"] - df["lag_p"]
df = df.dropna(subset=["jump"])
thresh = df["jump"].quantile(0.90)
tr = df[df.jump >= thresh].groupby("ISO3").first().reset_index()
tr = tr[["ISO3", "Year"]].rename(columns={"Year": "event_year"})
country = df[["ISO3", "Country"]].drop_duplicates()
out = tr.merge(country, on="ISO3")[["ISO3", "Country", "event_year"]].sort_values("Country")
out.to_csv("fraser_treated_country_year.csv", index=False)
print(out.to_string(index=False))
