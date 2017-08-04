# W. SEWLAL 1383337
module FileIOExtensionsCPP

# using JLD2

export checkDirectory, makeDirectory, checkFile, loadFile, fnameNew, filePathCheck
export pathFigures, pathSavedZ, pathSavedResults, pathSavedLogs

const today = string(Dates.today())
projectFiles = ["CPPParallelTempering.jl", "CPPParallelTempering.jl"]
# projectFiles = ["ParametersCPPMH.jl", "ParallelTemperingCPPMH.jl", "moduleCheck.jl", "PlottingCPPMH.jl"]


"""
checkDirectory(dir::AbstractString)

Throws an error if the directory 'dir' does not exist
"""
function checkDirectory(dir::AbstractString)
    !isdir(dir) && throw(ArgumentError(string("Directory ", dir, " not found")))
    return
end

"""
checkFile(fname::AbstractString, dir::AbstractString)

Throws an error if file 'dir/fname' does not exist
"""
function checkFile(fname::AbstractString, dir::AbstractString)
    !isfile(joinpath(dir, fname)) && throw(ArgumentError(string("File", fname, " not found")))
    return
end

"""
makeDirectory(dir::AbstractString)

Creates 'dir' if it does not exists. The function will attempt to create the parent directories of 'dir' recursively if they do not exist.
"""
function makeDirectory(dir::AbstractString)
    if (!isdir(dir))
        makeParentDirectory(dir)
        mkdir(dir)
    end
    return
end

"""
makeParentDirectory(dir::AbstractString)

Recursively checks and, if needed, creates parent directories of 'dir'.
"""
function makeParentDirectory(dir::AbstractString)
    parentDir = dirname(dir)
    if (!isdir(parentDir))
        warn(string("Directory ", parentDir, " does not exists. Creating parent directory."))
        makeParentDirectory(parentDir)
        mkdir(parentDir)
    end
    return
end

function getLastFilename(fnameBase::AbstractString, fnameExtension::AbstractString, dir::AbstractString)
    checkDirectory(dir)
    counter::Int64 = 1
    fname::AbstractString = string(fnameBase, "-", counter, fnameExtension)
    checkFile(fname, dir)

    while isfile(joinpath(dir, fname))
        lastfname = joinpath(dir, fname)
        counter += 1
        fname = string(fnameBase, "-", counter, fnameExtension)
    end
    return lastfname
end
    # fname::AbstractString = string(fnameBase, fnameExtension)

function getFilename(fnameBase::AbstractString, fnameExtension::AbstractString, dir::AbstractString, counter::Int64)
    checkDirectory(dir)
    fname::AbstractString = string(fnameBase, "-", counter, fnameExtension)
    checkFile(fname, dir)
    return joinpath(dir, fname)
end

function getFilename(fnameBase::AbstractString, fnameExtension::AbstractString, dir::AbstractString)
    return getLastFilename(fnameBase, fnameExtension, dir)
end

function getLastFilenameCounter(fnameBase::AbstractString, fnameExtension::AbstractString, dir::AbstractString)
    checkDirectory(dir)
    counter::Int64 = 1
    fname::AbstractString = string(fnameBase, "-", counter, fnameExtension)

    while isfile(joinpath(dir, fname))
        counter += 1
        fname = string(fnameBase, "-", counter, fnameExtension)
    end
    return (counter - 1)
end


function proposeFilename(fnameBase::AbstractString, fnameExtension::AbstractString, dir::AbstractString)
    checkDirectory(dir)
    fname::AbstractString = string(fnameBase, "-", 1, fnameExtension) # not needed
    !isfile(joinpath(dir, fname)) && return joinpath(dir, fname) # not needed
    warn(string("File ", fname, " already exists.")) # not needed
    counter = getLastFilenameCounter(fnameBase, fnameExtension, dir) + 1
    fname = string(fnameBase, "-", counter, fnameExtension)    
    warn(string("Creating ", fname))
    return joinpath(dir, fname)
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