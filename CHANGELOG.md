# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.3.0 - 2024-09-24
### Changed
- Update mrdetect module to version 2.0.0 
- Split detectSNVs task into three components

## 1.2.0 - 2024-09-24
### Added
- Resource specific files now exposed in map and required inputs reference and instrument model determine resources to be used
- Get HBC BAM list path from environment module `pwgs-hbc` instead of hard-coded value
### Changed
- Update HBC BAM files to version 2.0

## 1.1.1 -2024-08-01
### Fixed
- Fixed output file name issues. [GRD-809](https://jira.oicr.on.ca/browse/GRD-809)

## 1.1.0 - 2024-06-25
### Added
- Add vidarr labels to outputs (changes to medata only) [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) 

## 1.0.3 - 2023-01-16
### Added
- Added SNP.count and VAF files
### Changed
- Moved mrdetect algorithm scripts to bitbucket
- Separated filterVCF function from detectSNVs
- Changed filterAndDetect command to match code updates (v1.1.1)

## 1.0.2 - 2022-11-21
### Added
- Workflow finalized for production

## 1.0.0 - 2022-10-14
### Added
- A brand-new workflow.
