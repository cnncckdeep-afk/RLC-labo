# =============================================================
# RLC 串联电路：幅频特性曲线 + 谐振特性曲线
# 数据：表3-19-1-测量I-f谐振曲线数据表(2).xlsx
# 仪器：电阻箱 R、电容箱 C、电感箱 L
#   由 U_Ro、U_Rb 计算 I；由 I 与 L、C 标称值反算 U_C、U_L
# =============================================================

library(readxl)
library(ggplot2)

# 英文字体统一为 Times New Roman（坐标轴、图例、公式、注记）
if (.Platform$OS.type == "windows") {
  windowsFonts(TNR = windowsFont("Times New Roman"))
  font_en <- "TNR"
} else {
  font_en <- "Times New Roman"
}
update_geom_defaults("text",  list(family = font_en))
update_geom_defaults("label", list(family = font_en))

# ---------- 1. 读取数据 ----------
excel_file <- "表3-19-1-测量I-f谐振曲线数据表(3).xlsx"

raw <- read_excel(excel_file, sheet = 1, skip = 4, col_names = FALSE)
dat <- data.frame(
  f   = as.numeric(raw[[1]]),   # kHz
  URo = as.numeric(raw[[2]]),   # mV，R0 = 100 Ω 上的电压
  URb = as.numeric(raw[[4]])    # mV，R0' = 220 Ω 上的电压
)

# 元件标称值（R、C、L 箱）
L  <- 0.100          # H
C  <- 0.0500e-6      # F
R1 <- 100            # Ω
R2 <- 220            # Ω
U  <- 0.900          # V

# 回路电流（mA）
dat$I1 <- dat$URo / R1
dat$I2 <- dat$URb / R2

# 电阻上的实测电压（V）
dat$UR <- dat$URo / 1000

# 由实测电流与 L、C 标称值反算电压（左图用 R1 这一组）
omega <- 2 * pi * dat$f * 1000
dat$UL <- dat$I1 * 1e-3 * omega * L          # V
dat$UC <- dat$I1 * 1e-3 / (omega * C)        # V

# ---------- 2. 颜色（换回 No.021 蓝—白—红配色的两端） ----------
col_R100 <- "#104680"   # 深蓝（No.021 配色 1）
col_R220 <- "#B72230"   # 深红（No.021 配色 9）
col_ideal <- "#DC6D57"  # 珊瑚红（No.021 配色 8，理想曲线）
col_band_blue <- "#B6D7E8"  # 浅蓝（No.021 配色 4，蓝线误差带）
col_band_red  <- "#F6B293"  # 浅橙（No.021 配色 7，红线误差带）
col_UC   <- col_R100
col_UL   <- col_R220
col_UR   <- col_ideal    # 第三条曲线沿用 No.021 的第三个强调色（珊瑚红）

# ---------- 3. 统一主题：B5 图幅、白底、超细方格 ----------
# 由 A1（841 × 594 mm）等比缩至 B5（250 × 176 mm），所有尺寸乘 sc
sc <- 250 / 841
theme_rlc <- theme_bw(base_size = 48 * sc, base_family = font_en) +
  theme(
    text               = element_text(family = font_en),
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.border       = element_rect(colour = "black", linewidth = 1.0 * sc),
    axis.line          = element_line(colour = "black", linewidth = 2.4 * sc,
                                      lineend = "square"),
    panel.grid.major   = element_line(colour = "grey22", linewidth = 0.55 * sc),
    panel.grid.minor   = element_line(colour = "grey80", linewidth = 0.40 * sc),
    axis.ticks         = element_line(colour = "black", linewidth = 0.6 * sc),
    axis.ticks.length  = grid::unit(4 * sc, "mm"),
    axis.text          = element_text(colour = "black", size = 44 * sc),
    axis.title         = element_text(colour = "black", size = 54 * sc),
    legend.position      = c(0.99, 0.99),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = "white", colour = "black",
                                        linewidth = 0.6 * sc),
    legend.margin        = margin(5 * sc, 8 * sc, 5 * sc, 8 * sc),
    legend.key           = element_blank(),
    legend.key.width     = grid::unit(24 * sc, "mm"),
    legend.key.height    = grid::unit(14 * sc, "mm"),
    legend.spacing.y     = grid::unit(3 * sc, "mm"),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 44 * sc),
    plot.margin          = margin(20 * sc, 28 * sc, 18 * sc, 20 * sc)
  )

# 超细网格：主格 0.2 kHz，次格 0.02 kHz
minor_breaks_x <- as.vector(outer(seq(1.2, 3.0, 0.2), seq(0.02, 0.18, 0.02), "+"))

# 曲线自然延伸到坐标轴左右边界 1.2–3.2 kHz
x_full <- seq(1.2, 3.2, length.out = 700)

# 自然三次样条：用于左图选定的代表点
smooth_xy <- function(x, y) {
  sp <- spline(x, y, xout = x_full, method = "natural")
  data.frame(x = sp$x, y = sp$y)
}

# 谐振公式拟合：I(f)=Im/sqrt(1+Q^2(f/f0-f0/f)^2)，可自然外推到坐标轴边界
fit_res <- function(f, I, Q0) {
  nls(I ~ Im / sqrt(1 + Q^2 * (f / f0 - f0 / f)^2),
      start = list(Im = max(I), f0 = 2.25, Q = Q0),
      data = data.frame(f = f, I = I))
}

# 拟合两条谐振曲线（左图标注 f0、Q 也要用 m1 的结果）
m1 <- fit_res(dat$f, dat$I1, Q0 = 11)
m2 <- fit_res(dat$f, dat$I2, Q0 = 6)

# ---------- 4. 图一：幅频特性 U_C、U_L、U_R ~ f ----------
# 挑选能反映趋势的代表性频点（避开重复测量的近似频点）
sel_f <- c(1.4, 1.6, 1.8, 2.0, 2.1, 2.15, 2.2, 2.25,
           2.3, 2.35, 2.4, 2.5, 2.7, 2.9, 3.1)
sel <- dat[dat$f %in% sel_f, ]

amp_long <- rbind(
  data.frame(f = sel$f, U = sel$UC, series = "U_C"),
  data.frame(f = sel$f, U = sel$UL, series = "U_L"),
  data.frame(f = sel$f, U = sel$UR, series = "U_R")
)
amp_long$series <- factor(amp_long$series, levels = c("U_C", "U_L", "U_R"))

amp_smooth <- do.call(rbind, lapply(split(amp_long, amp_long$series), function(d) {
  out <- smooth_xy(d$f, d$U)
  out$series <- d$series[1]
  out
}))

minor_x_amp <- as.vector(outer(seq(2.0, 2.45, 0.05), seq(0.01, 0.04, 0.01), "+"))
minor_y_amp <- as.vector(outer(seq(0, 11, 1), seq(0.2, 0.8, 0.2), "+"))

# ---- 曲线取值与右支交点（让指标线在碰到曲线处截止） ----
curve_y <- function(d, x) as.numeric(approx(d$x, d$y, xout = x)$y)
curve_x_right <- function(d, level) {
  d  <- d[order(d$x), ]
  ip <- which.max(d$y)
  dd <- d[ip:nrow(d), ]
  as.numeric(approx(dd$y, dd$x, xout = level)$y)
}
UC_curve <- subset(amp_smooth, series == "U_C")
UL_curve <- subset(amp_smooth, series == "U_L")
UR_curve <- subset(amp_smooth, series == "U_R")

# ---- 左图标注值 ----
# f_C、f_L 取两条平滑曲线的峰值位置；f0 取电流谐振拟合值
fC <- with(UC_curve, x[which.max(y)])
fL <- with(UL_curve, x[which.max(y)])
f0_amp <- as.numeric(coef(m1)["f0"])
Q1  <- 2 * pi * f0_amp * 1000 * L / (R1 + 20)   # 计入电感箱内阻 r_L = 20 Ω
QU  <- Q1 * U

y_fC <- curve_y(UC_curve, fC)
y_f0 <- curve_y(UC_curve, f0_amp)
y_fL <- curve_y(UL_curve, fL)
x_QU <- curve_x_right(UC_curve, QU)

# 竖直指标线的下端落在 U_R 曲线上（指标线碰到曲线即止，不穿越 U_R）
y_R_C <- curve_y(UR_curve, fC)
y_R_0 <- curve_y(UR_curve, f0_amp)
y_R_L <- curve_y(UR_curve, fL)

p_amp <- ggplot(amp_smooth, aes(x = x, y = y, colour = series)) +
  # 水平指标 QU：从左边框画到 U_C 曲线右支
  annotate("segment", x = 2.0, xend = x_QU, y = QU, yend = QU,
           colour = "grey30", linetype = "dashed", linewidth = 1.7 * sc) +
  # 输入电压 U 水平参考线（手绘图中的 U 虚线，贯穿全宽）
  annotate("segment", x = 2.0, xend = 2.5, y = U, yend = U,
           colour = "grey30", linetype = "dashed", linewidth = 1.7 * sc) +
  # 竖直指标：从 U_R 曲线画到对应曲线（碰到曲线即止）
  annotate("segment", x = fC, xend = fC, y = y_R_C, yend = y_fC,
           colour = col_UC, linetype = "dashed", linewidth = 1.7 * sc) +
  annotate("segment", x = f0_amp, xend = f0_amp, y = y_R_0, yend = y_f0,
           colour = "grey20", linetype = "longdash", linewidth = 1.7 * sc) +
  annotate("segment", x = fL, xend = fL, y = y_R_L, yend = y_fL,
           colour = col_UL, linetype = "dotdash", linewidth = 1.7 * sc) +
  geom_line(linewidth = 2.7 * sc, lineend = "round") +
  geom_point(data = amp_long, aes(x = f, y = U, colour = series),
             size = 7.5 * sc, shape = 15, stroke = 0) +
  # 横轴指标 f_C、f0、f_L（引线引到图内底部，避开曲线）
  annotate("segment", x = fC, xend = 2.210, y = 6.05, yend = 6.26,
           colour = col_UC, linewidth = 1.10 * sc) +
  annotate("segment", x = f0_amp, xend = 2.235, y = 6.05, yend = 6.26,
           colour = "grey20", linewidth = 1.10 * sc) +
  annotate("segment", x = fL, xend = 2.290, y = 6.05, yend = 6.26,
           colour = col_UL, linewidth = 1.10 * sc) +
  annotate("label", x = 2.210, y = 6.30, label = "f[C]", parse = TRUE,
           size = 15 * sc, colour = col_UC, fill = "white", linewidth = 0) +
  annotate("label", x = 2.235, y = 6.30, label = "f[0]", parse = TRUE,
           size = 15 * sc, colour = "grey20", fill = "white", linewidth = 0) +
  annotate("label", x = 2.290, y = 6.30, label = "f[L]", parse = TRUE,
           size = 15 * sc, colour = col_UL, fill = "white", linewidth = 0) +
  # 水平指标 QU
  annotate("label", x = 2.02, y = QU + 0.10, label = "Q*U", parse = TRUE,
           size = 15 * sc, colour = "black", fill = "white", linewidth = 0,
           hjust = 0, vjust = 0) +
  # 输入电压 U 标签（左侧靠框）
  annotate("label", x = 2.02, y = U + 0.10, label = "U", parse = FALSE,
           size = 15 * sc, colour = "black",
           fill = "white", linewidth = 0, hjust = 0, vjust = 0) +
  # 曲线直接标注（对照手绘图）：U_C、U_L 置于峰顶上方，U_R 置于峰顶下方空区
  annotate("label", x = 2.205, y = 11.10, label = "U[C]", parse = TRUE,
           size = 14 * sc, colour = col_UC, fill = "white", linewidth = 0) +
  annotate("label", x = 2.272, y = 11.10, label = "U[L]", parse = TRUE,
           size = 14 * sc, colour = col_UL, fill = "white", linewidth = 0) +
  annotate("label", x = 2.238, y = 0.32, label = "U[R]", parse = TRUE,
           size = 14 * sc, colour = col_UR, fill = "white", linewidth = 0) +
  scale_colour_manual(
    values = c("U_C" = col_UC, "U_L" = col_UL, "U_R" = col_UR),
    labels = c(bquote(U[C]), bquote(U[L]), bquote(U[R]))
  ) +
  guides(colour = guide_legend(
    override.aes = list(shape = 22, size = 13 * sc,
                        fill = c(col_UC, col_UL, col_UR),
                        colour = NA, linewidth = 0)
  )) +
  scale_x_continuous(
    limits = c(2.0, 2.5),
    breaks = seq(2.0, 2.5, 0.05),
    labels = function(x) sprintf("%.2f", x),
    minor_breaks = minor_x_amp,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 12),
    breaks = seq(0, 12, 1),
    minor_breaks = minor_y_amp,
    expand = c(0, 0)
  ) +
  labs(x = "f (kHz)", y = "U (V)") +
  theme_rlc

ggsave("RLC_amplitude_Uc_UL.png", p_amp,
       width = 250, height = 176, units = "mm",
       dpi = 300, bg = "white")

# ---------- 5. 图二：谐振特性 I ~ f（两种 R + 理想曲线） ----------
# 绝对理想谐振曲线：完全由标称参数按串联谐振公式计算，不含任何拟合/测量偏差
#   I(f) = U / sqrt(R^2 + (ωL - 1/(ωC))^2)，取 R = R1 + r_L（电感箱内阻 20 Ω）
r_L <- 20
omega_full <- 2 * pi * x_full * 1000
I_ideal <- U / sqrt((R1 + r_L)^2 + (omega_full * L - 1 / (omega_full * C))^2) * 1000  # mA

cur_smooth <- rbind(
  data.frame(x = x_full, y = predict(m1, data.frame(f = x_full)), series = "R0_100"),
  data.frame(x = x_full, y = predict(m2, data.frame(f = x_full)), series = "R0_220"),
  data.frame(x = x_full, y = I_ideal, series = "Ideal")
)
cur_smooth$series <- factor(cur_smooth$series,
                            levels = c("R0_100", "R0_220", "Ideal",
                                       "Band_100", "Band_220"))

# 纵坐标测量误差带：以拟合曲线值为基准，上下各 5%（如 3.834 → 3.642–4.026）
band_pct <- 0.05
y_m1 <- predict(m1, data.frame(f = x_full))
y_m2 <- predict(m2, data.frame(f = x_full))
cur_band <- rbind(
  data.frame(x = x_full,
             ymin = y_m1 * (1 - band_pct),
             ymax = y_m1 * (1 + band_pct),
             series = "R0_100"),
  data.frame(x = x_full,
             ymin = y_m2 * (1 - band_pct),
             ymax = y_m2 * (1 + band_pct),
             series = "R0_220")
)
cur_band$series <- factor(cur_band$series, levels = c("R0_100", "R0_220"))

cur_points <- rbind(
  data.frame(f = dat$f, y = dat$I1, series = "R0_100"),
  data.frame(f = dat$f, y = dat$I2, series = "R0_220")
)
cur_points$series <- factor(cur_points$series, levels = c("R0_100", "R0_220"))

# 图例占位数据：仅用于让两条误差带以透明色块出现在图例中（面板中不可见）
band_key <- data.frame(
  f = 1.2, y = 8,
  series = factor(c("Band_100", "Band_220"),
                  levels = c("R0_100", "R0_220", "Ideal", "Band_100", "Band_220"))
)

minor_y_res <- as.vector(outer(seq(0, 7, 1), seq(0.2, 0.8, 0.2), "+"))

# ---- 右图标注值 ----
# 通频带边界取自 表3-19-3（f1、f2 由半功率点测得）
f1_1 <- 2.1573;  f2_1 <- 2.3483   # R = 120 Ω（100 Ω + r_L）
f1_2 <- 2.06879; f2_2 <- 2.45     # R = 240 Ω（220 Ω + r_L）
f0_res <- sqrt(f1_1 * f2_1)       # 谐振频率
Im1 <- as.numeric(coef(m1)["Im"]); Im2 <- as.numeric(coef(m2)["Im"])
half1 <- Im1 / sqrt(2); half2 <- Im2 / sqrt(2)

# 指标线端点：竖直方向到曲线，水平方向到右半支曲线
f0_m1 <- as.numeric(coef(m1)["f0"]); f0_m2 <- as.numeric(coef(m2)["f0"])
y_f1_1 <- as.numeric(predict(m1, data.frame(f = f1_1)))
y_f1_2 <- as.numeric(predict(m2, data.frame(f = f1_2)))
y_f0   <- as.numeric(predict(m1, data.frame(f = f0_res)))
y_f2_1 <- as.numeric(predict(m1, data.frame(f = f2_1)))
y_f2_2 <- as.numeric(predict(m2, data.frame(f = f2_2)))
right_root <- function(m, level) {
  f0m <- as.numeric(coef(m)["f0"])
  uniroot(function(x) as.numeric(predict(m, data.frame(f = x))) - level,
          c(f0m + 1e-6, 3.2))$root
}
x_half1 <- right_root(m1, half1)
x_half2 <- right_root(m2, half2)

p_res <- ggplot(cur_smooth, aes(x = x, y = y, colour = series)) +
  # 水平指标：从左框画到右半支曲线
  annotate("segment", x = 1.2, xend = f0_m1, y = Im1, yend = Im1,
           colour = "grey40", linetype = "dashed", linewidth = 1.7 * sc) +
  annotate("segment", x = 1.2, xend = x_half1, y = half1, yend = half1,
           colour = "grey40", linetype = "dashed", linewidth = 1.7 * sc) +
  annotate("segment", x = 1.2, xend = f0_m2, y = Im2, yend = Im2,
           colour = "grey40", linetype = "dashed", linewidth = 1.7 * sc) +
  annotate("segment", x = 1.2, xend = x_half2, y = half2, yend = half2,
           colour = "grey40", linetype = "dashed", linewidth = 1.7 * sc) +
  # 竖直指标：从下框画到曲线
  annotate("segment", x = f1_2, xend = f1_2, y = 0, yend = y_f1_2,
           colour = col_R220, linetype = "dashed", linewidth = 1.7 * sc) +
  annotate("segment", x = f1_1, xend = f1_1, y = 0, yend = y_f1_1,
           colour = col_R100, linetype = "dashed", linewidth = 1.7 * sc) +
  annotate("segment", x = f0_res, xend = f0_res, y = 0, yend = y_f0,
           colour = "grey20", linetype = "longdash", linewidth = 1.7 * sc) +
  annotate("segment", x = f2_1, xend = f2_1, y = 0, yend = y_f2_1,
           colour = col_R100, linetype = "dashed", linewidth = 1.7 * sc) +
  annotate("segment", x = f2_2, xend = f2_2, y = 0, yend = y_f2_2,
           colour = col_R220, linetype = "dashed", linewidth = 1.7 * sc) +
  # 误差带：上下限点连线并以浅实色填充，置于曲线下方
  geom_ribbon(data = cur_band,
              aes(x = x, ymin = ymin, ymax = ymax, fill = series),
              inherit.aes = FALSE,
              alpha = 0.35, colour = NA, show.legend = FALSE) +
  # key_glyph = "point" 让理想曲线（无数据点）也能在图例中显示色块
  geom_line(linewidth = 2.7 * sc, lineend = "round", key_glyph = "point") +
  geom_point(data = cur_points, aes(x = f, y = y, colour = series),
             size = 7.5 * sc, shape = 15, stroke = 0) +
  # 不可见的占位点，用于在图例中绘制透明误差带色块
  geom_point(data = band_key, aes(x = f, y = y, colour = series),
             size = 0, alpha = 0, stroke = 0) +
  # 横轴指标
  annotate("label", x = f1_2, y = 0.30, label = "f[1]*minute", parse = TRUE,
           size = 14.5 * sc, colour = col_R220, fill = "white", linewidth = 0) +
  annotate("label", x = f1_1, y = 0.30, label = "f[1]", parse = TRUE,
           size = 14.5 * sc, colour = col_R100, fill = "white", linewidth = 0) +
  annotate("label", x = f0_res, y = 0.30, label = "f[0]", parse = TRUE,
           size = 14.5 * sc, colour = "grey20", fill = "white", linewidth = 0) +
  annotate("label", x = f2_1, y = 0.30, label = "f[2]", parse = TRUE,
           size = 14.5 * sc, colour = col_R100, fill = "white", linewidth = 0) +
  annotate("label", x = f2_2, y = 0.30, label = "f[2]*minute", parse = TRUE,
           size = 14.5 * sc, colour = col_R220, fill = "white", linewidth = 0) +
  # 水平指标
  annotate("label", x = 1.235, y = Im1 + 0.12, label = "I[m]", parse = TRUE,
           size = 14 * sc, colour = "black", fill = "white", linewidth = 0,
           hjust = 0, vjust = 0) +
  annotate("label", x = 1.235, y = half1 + 0.12, label = "I[m]/sqrt(2)", parse = TRUE,
           size = 14 * sc, colour = "black", fill = "white", linewidth = 0,
           hjust = 0, vjust = 0) +
  annotate("label", x = 1.235, y = Im2 + 0.12, label = "I[m]*minute", parse = TRUE,
           size = 14 * sc, colour = "black", fill = "white", linewidth = 0,
           hjust = 0, vjust = 0) +
  annotate("label", x = 1.235, y = half2 + 0.12, label = "I[m]*minute/sqrt(2)",
           parse = TRUE, size = 14 * sc, colour = "black", fill = "white",
           linewidth = 0, hjust = 0, vjust = 0) +
  # 红色参数注记引线（仅一条，从红曲线引至红字块左缘）
  annotate("segment", x = 2.50, xend = 2.794,
           y = as.numeric(predict(m2, data.frame(f = 2.50))), yend = 4.18,
           colour = col_R220, linewidth = 1.1 * sc, lineend = "round") +
  # 拟合公式（横排，右对齐；左移 12 小格后再右移 6 小格）
  annotate("text", x = 2.917, y = 6.05,
           label = "I(f)==frac(I[m], sqrt(1+Q^2*(f/f[0]-f[0]/f)^2))",
           parse = TRUE, size = 12 * sc, colour = "black", hjust = 1) +
  # 拟合参数：R0、Im、f0、Q 四个指标各占一行；行距压到最小，分置蓝、红两框
  # 蓝色参数注记引线（从蓝曲线上升支引至蓝框右下角）
  annotate("segment", x = 2.05, xend = 1.95,
           y = as.numeric(predict(m1, data.frame(f = 2.05))), yend = 6.00,
           colour = col_R100, linewidth = 1.1 * sc, lineend = "round") +
  # 蓝字块（左上空白区）：白底黑框，先画矩形以遮住网格
  annotate("rect", xmin = 1.72, xmax = 1.95, ymin = 6.00, ymax = 7.00,
           fill = "white", colour = "black", linewidth = 0.6 * sc) +
  annotate("text", x = 1.746, y = 6.85, label = "R[0]==100~Omega",
           parse = TRUE, size = 12 * sc, colour = col_R100, hjust = 0) +
  annotate("text", x = 1.746, y = 6.62,
           label = sprintf("I[m]==%.3f~mA", coef(m1)["Im"]),
           parse = TRUE, size = 12 * sc, colour = col_R100, hjust = 0) +
  annotate("text", x = 1.746, y = 6.39,
           label = sprintf("f[0]==\"%.3f\"~kHz", coef(m1)["f0"]),
           parse = TRUE, size = 12 * sc, colour = col_R100, hjust = 0) +
  annotate("text", x = 1.746, y = 6.16,
           label = sprintf("Q==%.2f", coef(m1)["Q"]),
           parse = TRUE, size = 12 * sc, colour = col_R100, hjust = 0) +
  # 红字块（右侧红色块）：白底黑框
  annotate("rect", xmin = 2.794, xmax = 3.026, ymin = 3.454, ymax = 4.451,
           fill = "white", colour = "black", linewidth = 0.6 * sc) +
  annotate("text", x = 2.82, y = 4.30, label = "R[0]*minute==220~Omega",
           parse = TRUE, size = 12 * sc, colour = col_R220, hjust = 0) +
  annotate("text", x = 2.82, y = 4.07,
           label = sprintf("I[m]==%.3f~mA", coef(m2)["Im"]),
           parse = TRUE, size = 12 * sc, colour = col_R220, hjust = 0) +
  annotate("text", x = 2.82, y = 3.84,
           label = sprintf("f[0]==\"%.3f\"~kHz", coef(m2)["f0"]),
           parse = TRUE, size = 12 * sc, colour = col_R220, hjust = 0) +
  annotate("text", x = 2.82, y = 3.61,
           label = sprintf("Q==%.2f", coef(m2)["Q"]),
           parse = TRUE, size = 12 * sc, colour = col_R220, hjust = 0) +
  scale_fill_manual(values = c("R0_100" = col_band_blue,
                               "R0_220" = col_band_red)) +
  scale_colour_manual(
    values = c("R0_100" = col_R100, "R0_220" = col_R220, "Ideal" = col_ideal,
               "Band_100" = col_band_blue, "Band_220" = col_band_red),
    breaks = c("R0_100", "R0_220", "Ideal", "Band_100", "Band_220"),
    drop = FALSE,
    labels = c(bquote(R[0] == 100 ~ Omega ~ "(NLLS)"),
               bquote(R[0] * minute == 220 ~ Omega ~ "(NLLS)"),
               "Ideal curve (TC)",
               bquote('Error band (precision '*delta*' ±5%, '*R[0] == 100~Omega*')'),
               bquote('Error band (precision '*delta*' ±5%, '*R[0]*minute == 220~Omega*')'))
  ) +
  guides(colour = guide_legend(
    override.aes = list(shape = 22, size = 13 * sc, alpha = 1,
                        fill = c(col_R100, col_R220, col_ideal,
                                 adjustcolor(col_band_blue, alpha.f = 0.35),
                                 adjustcolor(col_band_red, alpha.f = 0.35)),
                        colour = NA, linewidth = 0)
  )) +
  scale_x_continuous(
    limits = c(1.2, 3.2),
    breaks = seq(1.2, 3.2, 0.2),
    minor_breaks = minor_breaks_x,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 8),
    breaks = seq(0, 8, 1),
    minor_breaks = minor_y_res,
    expand = c(0, 0)
  ) +
  labs(x = "f (kHz)", y = "I (mA)") +
  theme_rlc

ggsave("RLC_resonance_I_f.png", p_res,
       width = 250, height = 176, units = "mm",
       dpi = 300, bg = "white")

cat("完成：RLC_amplitude_Uc_UL.png, RLC_resonance_I_f.png\n")
