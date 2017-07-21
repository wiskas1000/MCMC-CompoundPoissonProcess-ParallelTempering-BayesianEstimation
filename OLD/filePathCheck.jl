# W. SEWLAL 1383337
module fileCheck

using JLD

export checkDirectory, makeDirectory, checkFile, loadFile, fnameNew, filePathCheck
export pathFigures, pathSavedZ, pathSavedResults, pathSavedLogs

const today = string(Dates.today())
projectFiles = ["ParametersCPPMH.jl", "ParallelTemperingCPPMH.jl", "moduleCheck.jl", "PlottingCPPMH.jl"]


function checkDirectory(dir::AbstractString)
    !isdir(dir) && throw(ArgumentError(string("Directory ", dir, " not found")))
    return
end


function makeDirectory(dir::AbstractString)
    !isdir(dir) && mkdir(dir)
    return
end


function checkFile(fname::AbstractString, dir::AbstractString)
    !isfile(joinpath(dir, fname)) && throw(ArgumentError(string("File", fname, " not found")))
    return
end

### Possible extension: fnameZ = latestFName(fnameZbase, fnameZextension, savedZPath), which 1. checkFile(fnameZ, savedZPath), 2. gives a warning if multiple savefiles are present and 3. returns the last savefile filename.
function loadFile(fnameBase::AbstractString, fnameExtension::AbstractString, dir::AbstractString)
    checkDirectory(dir)
    fname::AbstractString = string(fnameBase, fnameExtension)
    checkFile(fname, dir)
    return load(joinpath(dir, fname))
end


function fnameNew(fnameBase::AbstractString, fnameExtension::AbstractString, dir::AbstractString)
    checkDirectory(dir)
    fname::AbstractString = string(fnameBase, fnameExtension)
    !isfile(joinpath(dir, fname)) && return joinpath(dir, fname)
    warn(string("File ", fname, " already exists."))
    counter::Int64 = 1
    fname = string(fnameBase, "-", counter, fnameExtension)
    while isfile(joinpath(dir, fname))
        counter += 1
        fname = string(fnameBase, "-", counter, fnameExtension)
    end
    warn(string("Creating ", fname))
    return joinpath(dir, fname)
end


function filePathCheck(pathSavefiles::AbstractString, pathProject::AbstractString, savePlots::Bool)
    # Check savefile path
    !isdir(pathSavefiles) && throw(ArgumentError(string("pathSavefiles not set properly. Directory ", pathSavefiles, " not found")))
    global const pathFiguresBase = joinpath(pathSavefiles, "Figures")
    global const pathFigures = joinpath(pathFiguresBase, today)
    global const pathSavedZ = joinpath(pathSavefiles, "SimulationVariables")
    global const pathSavedResults = joinpath(pathSavefiles, "SimulationResults")
    global const pathSavedLogs = joinpath(pathSavefiles, "SimulationLogs")

    # Check/create savefile directories
    makeDirectory(pathFiguresBase)
    savePlots && makeDirectory(pathFigures)
    makeDirectory(pathSavedZ)
    makeDirectory(pathSavedResults)
    makeDirectory(pathSavedLogs)
    for dir in [pathSavefiles, pathFigures, pathSavedZ, pathSavedResults, pathSavedLogs]
        checkDirectory(dir)
    end

    # Check if projectfiles are present
    for file in projectFiles
        checkFile(file, pathProject)
    end
    # return
end


end
