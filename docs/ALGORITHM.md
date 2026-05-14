# Algorithm details

This document describes the math used by `gnss_nav` in enough detail that you
can (a) walk through it in class, and (b) extend it later. It mirrors the
flow described in `README.md`.

## 1. The pseudorange equation

The fundamental measurement we have is the **pseudorange** ρᵢ to satellite *i*:

$$
\rho_i = \|\mathbf{s}_i - \mathbf{r}\| + c(\delta t_r - \delta t_i^{sv}) + I_i + T_i + \varepsilon_i
$$

where

| Symbol | Meaning |
| --- | --- |
| **sᵢ** | SV ECEF position at *transmit* time |
| **r** | Receiver ECEF position |
| δt_r | Receiver clock offset (s) |
| δt_i^sv | SV clock offset (s) |
| I_i, T_i | Iono / tropo path delay (m) |
| ε_i | Multipath + receiver noise |

We treat δt_r as unknown (one global scalar per epoch); everything else we
either measure, model, or compute.

## 2. SV position from ephemeris (IS-GPS-200, §20.3.3.4.3)

Given the broadcast Keplerian elements, the procedure is:

1. **Time from ephemeris**: tk = t − t_oe (wrapped to ±302400 s).
2. **Mean motion**: n = √(μ/a³) + Δn, with a = (√a)².
3. **Mean anomaly**: M = M₀ + n·tk.
4. **Eccentric anomaly** E: solve M = E − e·sin E (Newton iteration).
5. **True anomaly** ν: tan(ν/2) = √((1+e)/(1−e)) · tan(E/2).
6. **Argument of latitude** φ = ν + ω, then add harmonic perturbations:
   - δu = Cus sin 2φ + Cuc cos 2φ
   - δr = Crs sin 2φ + Crc cos 2φ
   - δi = Cis sin 2φ + Cic cos 2φ
7. Updated u, r, i.
8. **Orbital plane position**: x' = r cos u, y' = r sin u.
9. **Corrected longitude of ascending node** including Earth rotation:
   Ω = Ω₀ + (Ω̇ − ω_e)·tk − ω_e·t_oe.
10. Rotate orbital → ECEF using Ω and i.

Velocities are the analytic derivatives of each of the above quantities — see
`compute_sv_state` for the full expansion.

## 3. SV clock correction

$$
\delta t^{sv}_i = a_{f0} + a_{f1}(t - t_{oc}) + a_{f2}(t - t_{oc})^2 + \Delta t_r
$$

where the **relativistic** correction is

$$
\Delta t_r = F \cdot e \cdot \sqrt{a} \cdot \sin E, \quad F = -4.442807633 \times 10^{-10}.
$$

We apply `+c·δt^sv` to the pseudorange (i.e. add it; the SV clock being fast
makes the pseudorange short).

## 4. Earth rotation

The SV is transmitting in an Earth-fixed frame that has rotated by ω_e·τ during
the signal's flight (τ ≈ 70 ms). We rotate the SV position by −ω_e·τ around Z
so it sits in the ECEF frame at *reception* time. This is folded into the
ephemeris computation via the corrected Ω above (term `−OMEGA_E_DOT * (toe + tau)`).

## 5. Atmospheric corrections

- **Klobuchar** (`corrections.klobuchar_delay`): single-frequency broadcast
  model using 8 coefficients from the nav file header. Returns L1 delay,
  scaled to the actual carrier by 1/f².
- **Saastamoinen** (`corrections.saastamoinen_delay`): tropospheric delay with
  a standard atmosphere model. No surface met required.

## 6. Least-squares position solver

After corrections the residual equation for each SV is

$$
\rho_i - \|\mathbf{s}_i - \mathbf{r}\| - c\delta t_r = \varepsilon_i
$$

We linearize around an initial guess **r₀**:

$$
\rho_i - \rho_i^{(0)} \approx -\hat{e}_i^T \Delta\mathbf{r} + c\Delta\delta t_r
$$

where ê_i is the unit line-of-sight from receiver to SV. Stacking over all
SVs gives **HΔx = δρ**, where

$$
H = \begin{bmatrix} -\hat{e}_1^T & 1 \\ \vdots & \vdots \\ -\hat{e}_n^T & 1 \end{bmatrix}.
$$

We solve the weighted normal equations HᵀWH·Δx = HᵀW·δρ and iterate
Gauss-Newton until ‖Δr‖ < 1 mm or max-iter is reached. Weights are sin²(elev).

**GDOP** is √trace((HᵀH)⁻¹).

## 7. Doppler-based velocity

For each SV the range-rate is

$$
\dot\rho_i = -\lambda_i D_i = (\mathbf{v}_i^{sv} - \mathbf{v}_r) \cdot \hat{e}_i + c\dot{\delta t}_r
$$

Stacked LSQ gives **(v_r, c·δṫ_r)** in one shot. We use the same weights as
the position solver. If Doppler isn't available we fall back to differencing
two consecutive positions divided by epoch dt.

## 8. Time and output

- RINEX records receiver time stamps that we treat as GPS-system time and
  also convert to UTC by subtracting the GPS-UTC leap offset
  (`GPS_UTC_LEAP_SECONDS = 18` as of 2017).
- ECEF → WGS-84 LLA uses Bowring iteration in `ecef_to_lla`.
- CSV: one row per fix. KML: a `<LineString>` for the path plus per-fix
  `<Placemark>` with `<TimeStamp>` so Google Earth can animate.
