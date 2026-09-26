"""
SED plotting module for single and binary stellar systems.
"""

import numpy as np
import matplotlib.pyplot as plt


def _extrapolate_spectrum(wave, flux, wave_grid):
    """
    Interpolate and extrapolate spectrum to a wavelength grid.

    Uses power-law extrapolation for wavelengths beyond the model range.
    """
    max_wave = wave.max()
    flux_grid = np.zeros_like(wave_grid, dtype=float)

    # Interpolate within range
    in_range = wave_grid <= max_wave
    flux_grid[in_range] = np.interp(wave_grid[in_range], wave, flux)

    # Extrapolate beyond range using power law fit to red end
    if (~in_range).any():
        red_mask = wave > 0.8 * max_wave
        if red_mask.sum() > 2:
            log_w = np.log10(wave[red_mask])
            log_f = np.log10(np.maximum(flux[red_mask], 1e-30))
            slope, intercept = np.polyfit(log_w, log_f, 1)
            flux_grid[~in_range] = 10**(slope * np.log10(wave_grid[~in_range]) + intercept)

    return flux_grid


def plot_single_sed(
    result,
    photometry_df,
    title=None,
    figsize=(14, 9),
    xlim=(1000, 1e5),
    ylim=None,
    residual_ylim=(-5, 5),
    show_filter_labels=False,
    show_residuals=True,
    fill_under=False,
    fill_alpha=0.3,
    source_label=None,
    ax=None,
    fig=None,
    data_kwargs=None,
    model_kwargs=None,
    fill_kwargs=None,
    residual_mode='sigma',
    residual_offset=0.0,
    residual_ylabel=None,
):
    """
    Plot a single star SED fit result.

    Parameters
    ----------
    result : dict
        Result dictionary from SingleSEDFitter.fit() or fit_single_sed().
    photometry_df : pandas.DataFrame
        Photometric data used for fitting.
    title : str, optional
        Plot title. If None, auto-generates from fit parameters.
    figsize : tuple, optional
        Figure size (default: (14, 9)).
    xlim : tuple, optional
        X-axis limits in Angstroms (default: (1000, 1e5)).
    ylim : tuple, optional
        Y-axis limits for flux. If None, auto-scales.
    residual_ylim : tuple, optional
        Y-axis limits for residual panel (default: (-5, 5)).
    show_filter_labels : bool, optional
        Annotate data points with filter names (default: False).
    show_residuals : bool, optional
        Show residual panel below main plot (default: True).
    fill_under : bool, optional
        Fill under the model curve (default: False).
    fill_alpha : float, optional
        Alpha for fill_between (default: 0.3).
    ax : matplotlib.axes.Axes or tuple, optional
        If provided, plot on this axes. For residuals, pass (ax_main, ax_resid).
    fig : matplotlib.figure.Figure, optional
        If provided with ax, use this figure.
    data_kwargs : dict, optional
        Additional kwargs passed to errorbar for observed data.
    model_kwargs : dict, optional
        Additional kwargs passed to plot for model spectrum.
    fill_kwargs : dict, optional
        Additional kwargs passed to fill_between for model fill.
    residual_mode : str, optional
        How to compute residuals. Options:
        - 'sigma' (default): (obs - model) / error, normalized by uncertainty
        - 'absolute': (obs - model) - offset, raw distance from offset
    residual_offset : float, optional
        When residual_mode='absolute', sets symmetric y-limits to (-offset, +offset).
        If 0.0, uses residual_ylim instead. Ignored when residual_mode='sigma'.
    residual_ylabel : str, optional
        Custom y-axis label for residuals. If None, uses default labels based on mode.

    Returns
    -------
    tuple
        (fig, axes) where axes is (ax_main, ax_resid) or ax_main if no residuals.
    """
    # Set default kwargs
    data_kwargs = data_kwargs or {}
    model_kwargs = model_kwargs or {}
    fill_kwargs = fill_kwargs or {}

    # Get spectrum and normalization
    spec = result['spectrum']
    norm = result['norm']
    model_wave = spec['wavelength'].values
    model_flux = spec['flux'].values * norm

    obs_wave = photometry_df['sed_wave'].values
    obs_flux = photometry_df['sed_flux'].values
    obs_eflux = photometry_df['sed_eflux'].values

    # Create figure if not provided
    if ax is None:
        if show_residuals:
            fig, (ax1, ax2) = plt.subplots(
                2, 1, figsize=figsize, height_ratios=[3, 1], sharex=True
            )
        else:
            fig, ax1 = plt.subplots(figsize=figsize)
            ax2 = None
    else:
        if isinstance(ax, tuple):
            ax1, ax2 = ax
        else:
            ax1 = ax
            ax2 = None
        if fig is None:
            fig = ax1.figure

    # Default data plotting options
    data_defaults = dict(
        fmt='o', markersize=8, label='Observed',
        capsize=3, color='black', zorder=10
    )
    data_defaults.update(data_kwargs)

    # Plot observed data
    ax1.errorbar(obs_wave, obs_flux, yerr=obs_eflux, **data_defaults)

    # Default model plotting options
    if source_label is None:
        source_label = f"Model (Teff={result['teff']}K, log g={result['logg']})"
    model_defaults = dict(
        color='r', alpha=0.7, linewidth=1,
        label=source_label
    )
    model_defaults.update(model_kwargs)

    # Plot model spectrum
    ax1.plot(model_wave, model_flux, **model_defaults)

    # Fill under model curve
    if fill_under:
        fill_defaults = dict(
            alpha=fill_alpha,
            color=model_defaults.get('color', 'r'),
        )
        fill_defaults.update(fill_kwargs)
        # For log scale, fill down to a small positive value
        ax1.fill_between(model_wave, model_flux, 1e-30, **fill_defaults)

    # Add filter labels if requested
    if show_filter_labels and 'sed_filter' in photometry_df.columns:
        for _, row in photometry_df.iterrows():
            ax1.annotate(
                row['sed_filter'].split(':')[-1],
                (row['sed_wave'], row['sed_flux']),
                textcoords='offset points', xytext=(5, 5),
                fontsize=8, alpha=0.7
            )

    ax1.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.legend()
    ax1.set_xlim(xlim)
    if ylim:
        ax1.set_ylim(ylim)

    # Auto-generate title if not provided
    if title is None:
        title = (
            f"Best-fit Model | Teff={result['teff']}K, log g={result['logg']}, "
            f"R={result['radius_rsun']:.4f} R_sun"
        )
    ax1.set_title(title)

    # Plot residuals
    if ax2 is not None:
        model_flux_interp = np.interp(obs_wave, model_wave, model_flux)

        if residual_mode == 'sigma':
            residuals = (obs_flux - model_flux_interp) / obs_eflux
            yerr = 1
            default_ylabel = r'Residuals ($\sigma$)'
            ylim = residual_ylim
        else:  # 'absolute'
            residuals = obs_flux - model_flux_interp
            yerr = obs_eflux
            default_ylabel = 'Residual'
            # Use residual_offset to set symmetric y-limits
            if residual_offset != 0.0:
                ylim = (-residual_offset, residual_offset)
            else:
                ylim = residual_ylim

        ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax2.errorbar(
            obs_wave, residuals, yerr=yerr,
            fmt='o', markersize=8, color='black', capsize=3
        )
        ax2.set_xlabel(r'Wavelength ($\mathring{\rm A}$)')
        ax2.set_ylabel(residual_ylabel if residual_ylabel is not None else default_ylabel)
        ax2.set_ylim(ylim)

        return fig, (ax1, ax2)

    return fig, ax1


def plot_binary_sed(
    result,
    photometry_df,
    title=None,
    figsize=(14, 9),
    xlim=(1000, 1e5),
    ylim=None,
    residual_ylim=(-5, 5),
    show_filter_labels=False,
    show_residuals=True,
    show_components=True,
    fill_under=False,
    fill_alpha=0.3,
    source1_label=None,
    source2_label=None,
    source1_color='blue',
    source2_color='orange',
    combined_color='red',
    ax=None,
    fig=None,
    data_kwargs=None,
    source1_kwargs=None,
    source2_kwargs=None,
    combined_kwargs=None,
    source1_fill_kwargs=None,
    source2_fill_kwargs=None,
    residual_mode='sigma',
    residual_offset=0.0,
    residual_ylabel=None,
):
    """
    Plot a binary SED fit result with individual components.

    Parameters
    ----------
    result : dict
        Result dictionary from BinarySEDFitter.fit() or fit_binary_sed().
    photometry_df : pandas.DataFrame
        Photometric data used for fitting.
    title : str, optional
        Plot title. If None, auto-generates from fit parameters.
    figsize : tuple, optional
        Figure size (default: (14, 9)).
    xlim : tuple, optional
        X-axis limits in Angstroms (default: (1000, 1e5)).
    ylim : tuple, optional
        Y-axis limits for flux. If None, auto-scales.
    residual_ylim : tuple, optional
        Y-axis limits for residual panel (default: (-5, 5)).
    show_filter_labels : bool, optional
        Annotate data points with filter names (default: False).
    show_residuals : bool, optional
        Show residual panel below main plot (default: True).
    show_components : bool, optional
        Plot individual source components (default: True).
    fill_under : bool, optional
        Fill under the component curves (default: False).
    fill_alpha : float, optional
        Alpha for fill_between (default: 0.3).
    source1_label : str, optional
        Legend label for source 1. If None, auto-generates.
    source2_label : str, optional
        Legend label for source 2. If None, auto-generates.
    source1_color : str, optional
        Color for source 1 component (default: 'blue').
    source2_color : str, optional
        Color for source 2 component (default: 'orange').
    combined_color : str, optional
        Color for combined model (default: 'red').
    ax : matplotlib.axes.Axes or tuple, optional
        If provided, plot on this axes. For residuals, pass (ax_main, ax_resid).
    fig : matplotlib.figure.Figure, optional
        If provided with ax, use this figure.
    data_kwargs : dict, optional
        Additional kwargs passed to errorbar for observed data.
    source1_kwargs : dict, optional
        Additional kwargs passed to plot for source 1 component.
    source2_kwargs : dict, optional
        Additional kwargs passed to plot for source 2 component.
    combined_kwargs : dict, optional
        Additional kwargs passed to plot for combined model.
    source1_fill_kwargs : dict, optional
        Additional kwargs passed to fill_between for source 1.
    source2_fill_kwargs : dict, optional
        Additional kwargs passed to fill_between for source 2.
    residual_mode : str, optional
        How to compute residuals. Options:
        - 'sigma' (default): (obs - model) / error, normalized by uncertainty
        - 'absolute': (obs - model) - offset, raw distance from offset
    residual_offset : float, optional
        When residual_mode='absolute', sets symmetric y-limits to (-offset, +offset).
        If 0.0, uses residual_ylim instead. Ignored when residual_mode='sigma'.
    residual_ylabel : str, optional
        Custom y-axis label for residuals. If None, uses default labels based on mode.

    Returns
    -------
    tuple
        (fig, axes) where axes is (ax_main, ax_resid) or ax_main if no residuals.
    """
    # Set default kwargs
    data_kwargs = data_kwargs or {}
    source1_kwargs = source1_kwargs or {}
    source2_kwargs = source2_kwargs or {}
    combined_kwargs = combined_kwargs or {}
    source1_fill_kwargs = source1_fill_kwargs or {}
    source2_fill_kwargs = source2_fill_kwargs or {}

    s1 = result['source1']
    s2 = result['source2']

    # Get spectra
    s1_spec = s1['spectrum']
    s2_spec = s2['spectrum']
    s1_norm = s1['norm']
    s2_norm = s2['norm']

    obs_wave = photometry_df['sed_wave'].values
    obs_flux = photometry_df['sed_flux'].values
    obs_eflux = photometry_df['sed_eflux'].values

    # Create wavelength grid for smooth plotting
    wave_grid = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 5000)

    # Compute flux on grid with extrapolation
    s1_wave = s1_spec['wavelength'].values
    s1_flux_raw = s1_spec['flux'].values
    s1_flux_grid = _extrapolate_spectrum(s1_wave, s1_flux_raw * s1_norm, wave_grid)

    s2_wave = s2_spec['wavelength'].values
    s2_flux_raw = s2_spec['flux'].values
    s2_flux_grid = _extrapolate_spectrum(s2_wave, s2_flux_raw * s2_norm, wave_grid)

    combined_flux_grid = s1_flux_grid + s2_flux_grid
    # Clip for log plot
    combined_flux_grid_plot = np.maximum(combined_flux_grid, 1e-30)

    # Create figure if not provided
    if ax is None:
        if show_residuals:
            fig, (ax1, ax2) = plt.subplots(
                2, 1, figsize=figsize, height_ratios=[3, 1], sharex=True
            )
        else:
            fig, ax1 = plt.subplots(figsize=figsize)
            ax2 = None
    else:
        if isinstance(ax, tuple):
            ax1, ax2 = ax
        else:
            ax1 = ax
            ax2 = None
        if fig is None:
            fig = ax1.figure

    # Default data plotting options
    data_defaults = dict(
        fmt='o', markersize=10, label='Observed',
        capsize=3, color='black', zorder=10
    )
    data_defaults.update(data_kwargs)

    # Plot observed data
    ax1.errorbar(obs_wave, obs_flux, yerr=obs_eflux, **data_defaults)

    # Plot individual components
    if show_components:
        # Source 1 defaults
        if source1_label is None:
            source1_label = f"Source 1 (Teff={s1['teff']}K, log g={s1['logg']})"
        s1_defaults = dict(alpha=0.6, linewidth=1.5, color=source1_color, label=source1_label)
        s1_defaults.update(source1_kwargs)
        ax1.plot(wave_grid, s1_flux_grid, '-', **s1_defaults)

        # Fill under source 1
        if fill_under:
            s1_fill_defaults = dict(
                alpha=fill_alpha,
                color=s1_defaults.get('color', source1_color),
            )
            s1_fill_defaults.update(source1_fill_kwargs)
            ax1.fill_between(wave_grid, s1_flux_grid, 1e-30, **s1_fill_defaults)

        # Source 2 defaults
        if source2_label is None:
            source2_label = f"Source 2 (Teff={s2['teff']}K, log g={s2['logg']})"
        s2_defaults = dict(alpha=0.6, linewidth=1.5, color=source2_color, label=source2_label)
        s2_defaults.update(source2_kwargs)
        s2_flux_plot = np.maximum(s2_flux_grid, 1e-30)
        ax1.plot(wave_grid, s2_flux_plot, '-', **s2_defaults)

        # Fill under source 2
        if fill_under:
            s2_fill_defaults = dict(
                alpha=fill_alpha,
                color=s2_defaults.get('color', source2_color),
            )
            s2_fill_defaults.update(source2_fill_kwargs)
            ax1.fill_between(wave_grid, s2_flux_plot, 1e-30, **s2_fill_defaults)

    # Combined model defaults
    combined_defaults = dict(alpha=0.8, linewidth=2, color=combined_color, label='Combined')
    combined_defaults.update(combined_kwargs)
    ax1.plot(wave_grid, combined_flux_grid_plot, '-', **combined_defaults)

    # Add filter labels if requested
    if show_filter_labels and 'sed_filter' in photometry_df.columns:
        for _, row in photometry_df.iterrows():
            ax1.annotate(
                row['sed_filter'].split(':')[-1],
                (row['sed_wave'], row['sed_flux']),
                textcoords='offset points', xytext=(5, 5),
                fontsize=8, alpha=0.7
            )

    ax1.set_ylabel(r'Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.legend()
    ax1.set_xlim(xlim)
    if ylim:
        ax1.set_ylim(ylim)

    # Auto-generate title if not provided
    if title is None:
        title = (
            f"Binary SED Fit | "
            f"S1: Teff={s1['teff']}K, R={s1['radius_rsun']:.3f} R_sun | "
            f"S2: Teff={s2['teff']}K, R={s2['radius_rsun']:.3f} R_sun"
        )
    ax1.set_title(title)

    # Plot residuals
    if ax2 is not None:
        # Interpolate combined flux to observed wavelengths
        combined_flux_interp = np.interp(obs_wave, wave_grid, combined_flux_grid)

        if residual_mode == 'sigma':
            residuals = (obs_flux - combined_flux_interp) / obs_eflux
            yerr = 1
            default_ylabel = r'Residuals ($\sigma$)'
            ylim = residual_ylim
        else:  # 'absolute'
            residuals = obs_flux - combined_flux_interp
            yerr = obs_eflux
            default_ylabel = 'Residual'
            # Use residual_offset to set symmetric y-limits
            if residual_offset != 0.0:
                ylim = (-residual_offset, residual_offset)
            else:
                ylim = residual_ylim

        ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax2.errorbar(
            obs_wave, residuals, yerr=yerr,
            fmt='o', markersize=8, color='black', capsize=3
        )
        ax2.set_xlabel(r'Wavelength ($\mathring{\rm A}$)')
        ax2.set_ylabel(residual_ylabel if residual_ylabel is not None else default_ylabel)
        ax2.set_ylim(ylim)

        return fig, (ax1, ax2)

    return fig, ax1


def plot_sed(result, photometry_df, **kwargs):
    """
    Convenience function that auto-detects single vs binary fit and plots.

    Parameters
    ----------
    result : dict
        Result dictionary from any fitter.
    photometry_df : pandas.DataFrame
        Photometric data used for fitting.
    **kwargs
        Additional arguments passed to plot_single_sed or plot_binary_sed.

    Returns
    -------
    tuple
        (fig, axes) where axes is (ax_main, ax_resid) or ax_main if no residuals.
    """
    # Detect binary vs single based on result structure
    if 'source1' in result and 'source2' in result:
        return plot_binary_sed(result, photometry_df, **kwargs)
    else:
        return plot_single_sed(result, photometry_df, **kwargs)
