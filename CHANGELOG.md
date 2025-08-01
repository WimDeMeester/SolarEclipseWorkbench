# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- Add LAST command to the scripts
- Execute external commands from scripts
- Installation using `pip install solareclipseworkbench`

### Changed
- (Describe any changes here)

### Fixed
- (Describe any bug fixes here)

---

## [1.1.1] - 2025-04-28

### Fixed

- Fix crash when reference moment not known

## [1.1.0] - 2025-03-25

Version 1.1.0 fixes some bugs, makes it possible to take pictures through a telescope and provides new scripting possibilities. Version of the partial solar eclipse of March 29, 2025.

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