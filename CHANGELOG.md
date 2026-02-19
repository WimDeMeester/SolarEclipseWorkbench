# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- **Location Search & Saved Locations in GUI**: The Location popup in the main GUI now includes the
  same saved-locations drop-down and address-search functionality.
  - Saved-locations drop-down populated from `~/.sew_wizard_config.json`; last-used location is
    automatically selected when the popup opens.
  - Address-search bar (requires `geopy`) geocodes any city, street, or landmark via Nominatim and
    fetches elevation from the Open-Elevation API in a background thread.
  - "Save Location" button to persist a named location for future sessions.
  - The map auto-updates whenever coordinates change (300 ms debounce); the manual "Plot" button has
    been removed.
- **Shared `location_ui` Module**: `ConfigManager`, `GeocodingWorker`, and the new `LocationWidget`
  are now maintained in a single `location_ui.py` module and reused by both the wizard (`wizard.py`)
  and the main GUI (`gui.py`), eliminating code duplication.

### Changed
- (Describe any changes here)

### Fixed
- (Describe any changes here)


## [1.4.0] - 2026-02-18

### Added
- **Script Generation Wizard (`sew_wizard`)**: New PyQt6-based 5-page wizard for automated eclipse photography script generation
  - Page 1: Eclipse configuration with date, location, and free geocoding service
  - Page 2: Camera equipment settings (name, ISO range, aperture, sync intervals)
  - Page 3: Phenomena selection (partial, diamond ring, Baily's beads, corona, prominences, chromosphere, earthshine, voice prompts, solar filter)
  - Page 4: Summary and script preview
- **Exposure Calculator**: Scientific exposure calculation based on Xavier Jubier's data
  - 2D interpolation using sun altitude (0-60°) and observer altitude (0-3000m)
  - 18 phenomenon exposure tables covering all eclipse phases
  - Support for ND 4.0 and ND 5.0 solar filters
  - Automatic ISO adjustment when exposures exceed hand-held limits (1/30s)
  - Realistic camera shutter speeds (standard 1/3-stop increments)
- **Comprehensive Partial Phase Coverage**: Generate 300+ shots automatically
  - All partial phase shots from C1→C2 and C3→C4
  - User-configurable intervals (time-based or magnitude-based)
  - Sun altitude filtering to skip shots when sun is below horizon
  - Applied to C1/C4 contact moments and all partial phases
- **Free Geocoding Service**: Convert addresses to coordinates without API keys
  - Nominatim (OpenStreetMap) for address → latitude/longitude
  - Open-Elevation API for altitude lookup
  - Background thread processing to prevent UI blocking
  - Visual feedback with status messages
- **Smart Totality Optimization**:
  - Adaptive corona shot intervals based on totality duration
  - Automatic gap filling between major phenomenon
  - Earthshine feasibility checking (only during totality with sufficient time)
  - Prominence shots in early totality
  - Chromosphere shots before C3
  - 10-second buffer zones to prevent command overlap
- **Camera-Specific Optimizations**:
  - Nikon burst mode: parameter = number of pictures (default: 30)
  - Canon burst mode: parameter = duration in seconds (default: 3)
  - Automatic detection based on camera name
- **Time Format Improvements**: Consistent h:mm:ss.0 format throughout all scripts
- **CSV Parsing**: Proper handling of commas in command descriptions for reliable import/export

### Changed
- **Script Parsing**: Updated `utils.py` and `scripts.py` to use Python csv module for robust parsing



## [1.3.0] - 2026-02-04

### Added
- Introduce `BaseCamera` and a `VirtualCamera` for simulator mode.
- Add `GPhotoCameraAdapter` and vendor adapters (`CanonCamera`, `NikonCamera`) to wrap gphoto2 cameras.
- `get_camera_dict(..., is_simulator=True)` returns a `VirtualCamera` for easy testing and demos.
- Add gphoto-style stubs on `VirtualCamera` (`get_config`, `set_config`, `get_storageinfo`, `exit`) for compatibility.
- Background probing of cameras off the GUI thread; model updates scheduled on the main thread to avoid UI freezes.
- Added defensive helpers and fallbacks for gphoto2 errors (storage/time/config) and a one-time reinitialisation retry.
- Low-level gphoto fallbacks for capture operations to improve reliability with real cameras.
- Tests and example script for the virtual camera added (`tests/test_virtual_camera.py`, example_scripts/testVirtualCamera.txt).
- Documented simulator CLI flag in README; added runtime prints/logging to aid diagnostics.
 - Use Astropy/IERS to compute Delta T (TT − UT1) for eclipse reference times
 - Keep CSV-based Delta T as a fallback if no internet is available; added robust parsing of `td_ge`/`t0` and safeguards when ephemeris files are unavailable.


## [1.2.4] - 2026-01-16

### Added
- Add LAST command to the scripts
- Execute external commands from scripts
- Installation using `pip install solareclipseworkbench`
- Calculate the reference moments in Solar Eclipse Workbench, not using an external library anymore.  The Besselian elements are taken from the Five Millennium Canon of Solar Eclipses, which is available at https://eclipse.gsfc.nasa.gov/SEcat5/SEcat5.html.
- Extra information about the eclipses (maximum duration, type of eclipse, etc.) is now available in the drop-down menu to select eclipses.
- Added first unittests for the solar_eclipse module.
- Calculate the Besselian elements directly in the code.
- Adapt solar radius to the most recent value of 959.95 ±0.05 arcseconds as found by the Besselian elements team (https://www.besselianelements.com/).


## [1.1.1] - 2025-04-28

### Fixed

- Fix crash when reference moment not known

## [1.1.0] - 2025-03-25

Version 1.1.0 fixes some bugs, makes it possible to take pictures through a telescope and provides new scripting possibilities. Version for the partial solar eclipse of March 29, 2025.

### Added
- Take pictures when a camera is attached to a telescope

### Fixed
- No longer counting down after the eclipse
- Fix astronomy import
- Scripts updates

## [1.0.0] - 2024-03-29

### Added
First version of Solar Eclipse Workbench. To be tested during the total solar eclipse of April 8, 2024 in Mexico, USA and Canada.

- Add logo
- Add poetry and installation instructions
- Placeholders for basic commands
- Added all needed camera methods
- Camera overview
- Fixes for the camera class
- Placeholder for location-related functions
- Fix get_camera_overview after testing
- Scheduling commands
- Calculate eclipse reference moments
- Documentation updates
- Extra calulations for the reference times
- Notifications
- Fix timezone
- Placeholder for GUI
- Handling of different eclipse types
- Add support for annular eclipses
- Script to convert Solar Eclipse Maestro scripts
- Countdown to reference moments
- Schedule tasks
- Documentation + corrected formatting of countdown clocks
- First version of camera widget
- Schedule tasks
- Update camera overview
- Documentation + clean-up
- Synchronisation of the cameras
- General script and documentation updates
- Simulator mode + job scheduling
- Add take_burst method
- Start simulation at given time relative to reference moment.
- Add camera_name to CameraSettings
- Start of simulation + Display scheduled jobs
- Added simulator icon
- Fix crash if camera is not known
- Visualisation of scheduled jobs
- Camera updates
- Apply time format from setting to scheduled jobs
- Show duration for total eclipses
- Proper alignment of scheduled jobs table cells
- Improved robustness + proper formatting
- Camera overview + Time formatting
- Add possibility to start up gui.py using parameters for location and eclipse date.
- Improved logging + Added setters for controller
- Saving & loading settings
- Extract camera overview as dictionary
- Save settings via toolbar icon
- Made file chooser more robust
- Fix cameras_sync command
- Documented UI functionality
- Calculate next 20 solar eclipses and display in drop-down menu
- Fix cancel button of SettingsPopup
- Disable/enable camera icon in UI toolbar
- Add output from logger to file
- Disconnect camera(s) when UI is closed
- Convert TAKEBST and TAKEBKT to Solar Eclipse Workbench command in convert_sem_files script.
- Scroll to jobs that are up next
- Add pkg-config to installation instruction for apple silicon