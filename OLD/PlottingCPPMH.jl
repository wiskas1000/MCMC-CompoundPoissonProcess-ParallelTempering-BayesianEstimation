# W. SEWLAL 1383337 # 2015-12-03
module PlottingCPPMH
using fileCheck
using ParametersCPPMH
using Distributions
using Cairo
using Gadfly
using StatsBase

export pdfPlotResults, writecsvMH

const colors = [colorant"black", colorant"purple", colorant"blue", colorant"red", colorant"green", colorant"orange", colorant"brown", colorant"red", colorant"black", colorant"purple", colorant"blue", colorant"red", colorant"green", colorant"orange", colorant"brown", colorant"red"]
const colornames = ["black", "purple", "blue", "red", "green", "orange", "brown", "red", "black", "purple", "blue", "red", "green", "orange", "brown", "red"]


function plotLegend(G, colors, colornames, labels)
        title = string("Legend")
        myplot = plot(y = [0],
                        Theme(default_color = colors[1], guide_title_position = :left, key_position = :left, key_label_font_size = 24pt),
                        Guide.manual_color_key("Legend", labels, colornames[1:G]),
                    );
end


function plotLegend(G, colors, labels)
        title = string("Legend")
        myplot = plot(y = [0],
                        Theme(default_color = colors[1], guide_title_position = :left, key_position = :left, key_label_font_size = 24pt),
                        Guide.manual_color_key("Gamma", labels[1:4], colornames[1:4]),
#                         CHANGE THIS, REVERT TO labels, COLORNAMES[1:G]
                    );
end


function correctOrder!(MH_mu, MH_psi)
    c = 1:size(MH_mu, 2)
    rows = [sub(MH_mu, i, c) for i = 1:size(MH_mu, 1)]
    p = sortperm(rows; by=mean)
    MH_mu[:, :] = MH_mu[p, :]
    MH_psi[:, :] = MH_psi[p, :]
    return
end

function correctOrder!(PT_mu, PT_psi, J, G, M_PT)
    tmpArray::Array{Float64,2} = zeros(J, G)
    for g in 1:G
        tmpArray = reshape(PT_mu[:, g, :], J, M_PT)
        c = 1:size(tmpArray, 2)
        rows = [sub(tmpArray, i, c) for i = 1:size(tmpArray, 1)]
        p = sortperm(rows; by=mean)

        PT_mu[:, g, :] = PT_mu[p, g, :]
        PT_psi[:, g, :] = PT_psi[p, g, :]
    end
    return
end


function extractData!(MH_labda, MH_psi, MH_mu, MH_tau, J, M_MH, MHparamArray)
    for m in 1:M_MH
        MH_tau[m] = MHparamArray[m].tau
        MH_labda[m] = MHparamArray[m].labda
        for j in 1:J
            MH_psi[j, m] = MHparamArray[m].psi[j]
            MH_mu[j, m] = MHparamArray[m].mu[j]
        end
    end
    correctOrder!(MH_mu, MH_psi)

    MH_labda
    MH_psi
    MH_mu
    MH_tau
    return
end

function extractData!(PT_labda, PT_psi, PT_mu, PT_tau, J, G, M_PT, PTparamArray)
    for g in 1:G
        for m in 1:M_PT
            PT_tau[g, m] = PTparamArray[g, m].tau
            PT_labda[g, m] = PTparamArray[g, m].labda
            for j in 1:J
                PT_psi[j, g, m] = PTparamArray[g, m].psi[j]
                PT_mu[j, g, m] = PTparamArray[g, m].mu[j]
            end
        end
    end

    correctOrder!(PT_mu, PT_psi, J, G, M_PT)
    PT_labda
    PT_psi
    PT_mu
    PT_tau
    return
end

function sortRealParamArray!(real_param)
    v = sortperm(real_param.mu)
    real_param.mu = real_param.mu[v]
    real_param.psi = real_param.psi[v]
    real_param.ksi = real_param.ksi[v]
    return
end

function setLabels(gamma, G)
    labels = Array(AbstractString, G)
    for g in 1:G
        labels[g] = string(gamma[g])
    end
    return labels
end


function setAxisMinMax(real_param, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau)
    axis_mm_labda::Array{Float64,1} = zeros(Float64, 2)
    axis_mm_psi::Array{Float64,1} = zeros(Float64, 2)
    axis_mm_mu::Array{Float64,1} = zeros(Float64, 2)
    axis_mm_tau::Array{Float64,1} = zeros(Float64, 2)
    axis_mm_acceptance::Array{Float64,1} = zeros(Float64, 2)

    axis_mm_acceptance[2] = 1

#     if real_param.labda == 3
#         axis_mm_labda = [0, 4]
#         axis_mm_psi = [0, 6]
#         axis_mm_mu = [-3, 3]
#         axis_mm_tau = [0, 2]
#     elseif real_param.labda == 6
#         axis_mm_labda = [4, 8]
#         axis_mm_psi = [0, 8]
#         axis_mm_mu = [-3, 3]
#         axis_mm_tau = [0, 2]
#     else
        axis_mm_labda = [min(0, min(minimum(MH_labda), minimum(PT_labda))), max(0, max(maximum(MH_labda), maximum(PT_labda)))]
        axis_mm_psi = [min(0, min(minimum(MH_psi), minimum(PT_psi))), max(0, max(maximum(MH_psi), maximum(PT_psi)))]
        axis_mm_mu = [min(0, min(minimum(MH_mu), minimum(PT_mu))), max(0, max(maximum(MH_mu), maximum(PT_mu)))]
        axis_mm_tau = [min(0, min(minimum(MH_tau), minimum(PT_tau))), max(0, max(maximum(MH_tau), maximum(PT_tau)))]
#     end
    return axis_mm_labda, axis_mm_psi, axis_mm_mu, axis_mm_tau, axis_mm_acceptance
end


function plotParameter(inputParameter, inputName, real_value, axis_mm, xlabel, M, thinning, showTrueValue::Bool = true)
    if thinning
        steps::Int64 = 10
    else
        steps = 1
    end
    # title = string("Values for ", inputName)
    myplot = plot(y = reshape(inputParameter, M), x = collect(steps:steps:steps*M),
                                    Geom.line,
                                    Theme(default_color = colors[1]),
                                    Guide.xlabel(xlabel),
                                    Guide.ylabel(inputName),
                                    # Scale.x_continuous(format=:plain),
                                    Coord.Cartesian(ymin = axis_mm[1], ymax = axis_mm[2]),
                                    );
    showTrueValue && (append!(myplot.layers, layer(yintercept = [real_value], Geom.hline(color = colorant"orange", size = 2mm))))
    return myplot
end

function plotParameters(inputParameters, inputName, real_values, axis_mm, xlabel, J, M, thinning)
    plots = Array(Gadfly.Plot, J)
    showTrueValues::Bool = true
    if J != length(real_values)
        showTrueValues = false
        real_values = zeros(Float64, J)
    end
    for j in 1:J
        plots[j] = plotParameter(inputParameters[j, :], string(inputName, j), real_values[j], axis_mm, xlabel, M, thinning, showTrueValues)
    end
    return plots
end


function plotParameterLayered(inputParameter, inputName, real_value, axis_mm, xlabel, G, M, colors, thinning, showTrueValue::Bool = true)
    if thinning
        steps::Int64 = 10
    else
        steps = 1
    end
#     myplot = plotParameter(inputParameter[1, :], inputName, real_value, axis_mm, xlabel, M, thinning, false)
    myplot = plotParameter(inputParameter[1, :], inputName, 0, axis_mm, xlabel, M, thinning, false)
    for g in 2:G
            append!(myplot.layers, layer(y = inputParameter[g, :], x = collect(steps:steps:steps*M), Geom.line, Theme(default_color = colors[g])))
    end
    showTrueValue && (append!(myplot.layers, layer(yintercept = [real_value], Geom.hline(color = colorant"orange", size = 2mm))))
    return myplot
end

function plotParametersLayered(inputParameters, inputName, real_values, axis_mm, xlabel, J, G, M, colors, thinning)
    plots = Array(Gadfly.Plot, J)
    showTrueValues::Bool = true
    if J != length(real_values)
        showTrueValues = false
        real_values = zeros(Float64, J)
    end

    for j in 1:J
        plots[j] = plotParameterLayered(reshape(inputParameters[j, :, :], G, M), string(inputName, j), real_values[j], axis_mm, xlabel, G, M, colors, thinning, showTrueValues)
    end
    return plots
end



function fillPlotArray(plot1, plot2J, plot3, plot4J, plot5, J)
    myplots::Array{Gadfly.Plot, 1} = Array(Gadfly.Plot, 0)
    push!(myplots, plot1)
    for j in 1:J
        push!(myplots, plot2J[j])
    end
    push!(myplots, plot3)
    for j in 1:J
        push!(myplots, plot4J[j])
    end
    push!(myplots, plot5)
    return myplots
end

function mergeESS(ESS1, ESS2J, ESS3, ESS4J, J)
#     ESS = Array(Array{Float64,1}, 3 + 2 * J)
    ESS::Array{Array{Float64,1}, 1} = Array(Array{Float64,1}, 0)
    push!(ESS, ESS1)
    for j in 1:J
        push!(ESS, ESS2J[j])
    end
    push!(ESS, ESS3)
    for j in 1:J
        push!(ESS, ESS4J[j])
    end
    return ESS
end



function getGamma1(input::Array{Float64,2}, M)
    return reshape(input[1, :], M)
end

function getGamma1(input::Array{Float64,3}, J, M)
    return reshape(input[:, 1, :], J, M)
end

function plotDensityMHPT(inputParameterMH::Array{Float64,1}, inputParameterPT::Array{Float64,1}, M_MH::Int64, M_PT::Int64, real_value::Float64, xlabel::AbstractString, burnin::Int64, burninlabel, colors, showTrueValue::Bool = true)
    inputParameterMH = inputParameterMH[(burnin + 1):M_MH]
    inputParameterPT = inputParameterPT[(burnin + 1):M_PT]
    myplot = plot(layer(x = inputParameterMH, Geom.density, Theme(default_color = colors[1])),
                    layer(x = inputParameterPT, Geom.density, Theme(default_color = colors[2])),
#                     Guide.xlabel(xlabel),
                    Guide.xlabel(string(xlabel, " (burnin = ", burninlabel, ")")),
                )
    showTrueValue && (append!(myplot.layers, layer(xintercept = [real_value], Geom.vline(color = colorant"orange"))))
    return myplot
end

function plotDensityMHPT(inputParameterMH::Array{Float64,2}, inputParameterPT::Array{Float64,2}, J::Int64, M_MH::Int64, M_PT::Int64, real_values, xlabel::AbstractString, burnin::Int64, burninlabel, colors)
    plots = Array(Gadfly.Plot, J)
    showTrueValues::Bool = true
    if J != length(real_values)
        showTrueValues = false
        real_values = zeros(Float64, J)
    end

    for j in 1:J
        plots[j] = plotDensityMHPT(reshape(inputParameterMH[j, :], M_MH), reshape(inputParameterPT[j, :], M_PT), M_MH, M_PT, real_values[j], string(xlabel, j), burnin, burninlabel, colors, showTrueValues)
    end
    return plots
end

function plotDensitiesMHPT(real_param, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, PT_acceptance, J, G, M_MH, M_PT, burnin, burninlabel, labels, colors, colornamesMHPT, thinning)
    plot1 = plotDensityMHPT(MH_tau, getGamma1(PT_tau, M_PT), M_MH, M_PT, real_param.tau, "tau", burnin, burninlabel, colors)
    plot2J = plotDensityMHPT(MH_mu, getGamma1(PT_mu, J, M_PT), J, M_MH, M_PT, real_param.mu, "mu", burnin, burninlabel, colors)
    plot3 = plotDensityMHPT(MH_labda, getGamma1(PT_labda, M_PT), M_MH, M_PT, real_param.labda, "labda", burnin, burninlabel, colors)
    plot4J = plotDensityMHPT(MH_psi, getGamma1(PT_psi, J, M_PT), J, M_MH, M_PT, real_param.psi, "psi", burnin, burninlabel, colors)
    plot5 = plotLegend(2, colors, colornamesMHPT, ["MH (gamma = 1)", "PT (gamma = 1)"])
    return fillPlotArray(plot1, plot2J, plot3, plot4J, plot5, J)
end


# 
# FOR NOSWAP
# 
function plotDensityLayered(inputParameter, M::Int64, real_value::Float64, xlabel::AbstractString, burnin::Int64, burninlabel, colors, showTrueValue::Bool = true)
    inputParameter = inputParameter[:, (burnin + 1):M]
    myplot = plot(layer(x = inputParameter[1, :], Geom.density, Theme(default_color = colors[1])),
                    Guide.xlabel(string(xlabel, " (burnin = ", burninlabel, ")")),
                )
    G = length(inputParameter[:, 1])
    for g in 2:G
#             append!(myplot.layers, layer(x = inputParameter[g, :], Geom.density))
            append!(myplot.layers, layer(x = inputParameter[g, :], Geom.density, Theme(default_color = colors[g])))
    end
    showTrueValue && (append!(myplot.layers, layer(xintercept = [real_value], Geom.vline(color = colorant"orange"))))
    return myplot
end

function plotDensityLayered(inputParameter, J::Int64, M::Int64, real_values, xlabel::AbstractString, burnin::Int64, burninlabel, colors)
    plots = Array(Gadfly.Plot, J)
    showTrueValues::Bool = true
    if J != length(real_values)
        showTrueValues = false
        real_values = zeros(Float64, J)
    end
    G = length(inputParameter[1, :, 1])

    for j in 1:J
        plots[j] = plotDensityLayered(reshape(inputParameter[j, :, :], (G, M)), M, real_values[j], string(xlabel, j), burnin, burninlabel, colors, showTrueValues)
    end
    return plots
end

function plotDensitiesLayered(real_param, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, PT_acceptance, J, G, M_MH, M, burnin, burninlabel, labels, colors, colornamesMHPT, thinning)
    plot1 = plotDensityLayered(PT_tau, M, real_param.tau, "tau", burnin, burninlabel, colors)
    plot2J = plotDensityLayered(PT_mu, J, M, real_param.mu, "mu", burnin, burninlabel, colors)
    plot3 = plotDensityLayered(PT_labda, M, real_param.labda, "labda", burnin, burninlabel, colors)
    plot4J = plotDensityLayered(PT_psi, J, M, real_param.psi, "psi", burnin, burninlabel, colors)
#     plot5 = plotLegend(2, colors, colornamesMHPT, ["MH (gamma = 1)", "PT (gamma = 1)"])
    plot5 = plotLegend(G, colors, labels)
    return fillPlotArray(plot1, plot2J, plot3, plot4J, plot5, J)
end




function plotAutocorrelationMHPT(inputParameterMH::Array{Float64,1}, inputParameterPT::Array{Float64,1}, M_MH::Int64, M_PT::Int64, burnin::Int64, burninlabel, colors)
    axis_mm = [0, 1]
    ac_MH = autocor(inputParameterMH[(burnin + 1):M_MH])
    ac_PT = autocor(inputParameterPT[(burnin + 1):M_PT])
    ESS_MH = length((burnin + 1):M_MH) / (1 + 2 * sum(ac_MH))
    ESS_PT = length((burnin + 1):M_PT) / (1 + 2 * sum(ac_PT))
#     ESS_MH = length((burnin + 1):M_MH) / (1 + 2 * sum(abs(ac_MH)))
#     ESS_PT = length((burnin + 1):M_PT) / (1 + 2 * sum(abs(ac_PT)))

    myplot = plot(layer(y = ac_MH, x = collect(1:length(ac_MH)), Geom.point, Theme(default_color = colors[1])),
                    layer(y = ac_PT, x = collect(1:length(ac_PT)), Geom.point, Theme(default_color = colors[2])),
                    Guide.xlabel(string("Lag (burnin = ", burninlabel, ")")), # Guide.xlabel("Lag"),
                    Guide.ylabel("acf"),
                    Coord.Cartesian(ymin = min(0, minimum(ac_PT)), ymax = 1),
                    )
#     return myplot
    return myplot, [ESS_MH, ESS_PT]
end

function plotAutocorrelationMHPT(inputParameterMH::Array{Float64,2}, inputParameterPT::Array{Float64,2}, J::Int64, M_MH::Int64, M_PT::Int64, burnin::Int64, burninlabel, colors)
    plots = Array(Gadfly.Plot, J)
    ESS = Array(Array{Float64,1}, J)
    for j in 1:J
        plots[j], ESS[j] = plotAutocorrelationMHPT(reshape(inputParameterMH[j, :], M_MH), reshape(inputParameterPT[j, :], M_PT), M_MH, M_PT, burnin, burninlabel, colors)
    end
    return plots, ESS
end


function plotAutocorrelationsPTMH(MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, PT_acceptance, J, G, M_MH, M_PT, burnin, burninlabel, labels, colors, thinning)
    plot1, ESS1 = plotAutocorrelationMHPT(MH_tau, getGamma1(PT_tau, M_PT), M_MH, M_PT, burnin, burninlabel, colors)
    plot2J, ESS2J = plotAutocorrelationMHPT(MH_mu, getGamma1(PT_mu, J, M_PT), J, M_MH, M_PT, burnin, burninlabel, colors)
    plot3, ESS3 = plotAutocorrelationMHPT(MH_labda, getGamma1(PT_labda, M_PT), M_MH, M_PT, burnin, burninlabel, colors)
    plot4J, ESS4J = plotAutocorrelationMHPT(MH_psi, getGamma1(PT_psi, J, M_PT), J, M_MH, M_PT, burnin, burninlabel, colors)
    plot5 = plotLegend(G, colors, labels)
    return fillPlotArray(plot1, plot2J, plot3, plot4J, plot5, J), mergeESS(ESS1, ESS2J, ESS3, ESS4J, J)
end
#     ESSLABELS = ["tau", "mu", "labda", "psi"]


function plotParametersPTAllGamma(real_param, PT_labda, PT_psi, PT_mu, PT_tau, PTacceptance, J, G, M_PT, burnin, axis_mm_labda, axis_mm_psi, axis_mm_mu, axis_mm_tau, axis_mm_acceptance, labels, colors, UpdateModusMH, noSwap, thinning)
    if noSwap
        xlabel = string("Solo MCMC chains Iteration")
    else
        xlabel = string(UpdateModusMH, " Iteration")
    end
    #     xlabel = string("PT Iteration") #     xlabel = string("PT Iteration (", UpdateModusMH, ")")
    plot1 = plotParameterLayered(PT_tau, "tau", real_param.tau, axis_mm_tau, xlabel, G, M_PT, colors, thinning)
    plot2J = plotParametersLayered(PT_mu, "mu", real_param.mu, axis_mm_mu, xlabel, J, G, M_PT, colors, thinning)
    plot3 = plotParameterLayered(PT_labda, "labda", real_param.labda, axis_mm_labda, xlabel, G, M_PT, colors, thinning)
    plot4J = plotParametersLayered(PT_psi, "psi", real_param.psi, axis_mm_psi, xlabel, J, G, M_PT, colors, thinning)
    plot5 = plotParameterLayered(PTacceptance, "acceptance", 0, axis_mm_acceptance, xlabel, G, length(PTacceptance[1, :]), colors, thinning, false)
    return fillPlotArray(plot1, plot2J, plot3, plot4J, plot5, J)
end


function plotParametersPTGamma1(real_param, PT_labda, PT_psi, PT_mu, PT_tau, PT_acceptance, J, M_PT, burnin, axis_mm_labda, axis_mm_psi, axis_mm_mu, axis_mm_tau, axis_mm_acceptance, labels, UpdateModusMH, thinning)
    xlabel = string(UpdateModusMH, "(gamma=1) Iteration")
    #     xlabel = string("PT Iteration (", UpdateModusMH, ")")
    M_acceptance = length(PT_acceptance[1, :])
    plot1 = plotParameter(getGamma1(PT_tau, M_PT), "tau", real_param.tau, axis_mm_tau, xlabel, M_PT, thinning)
    plot2J = plotParameters(getGamma1(PT_mu, J, M_PT), "mu", real_param.mu, axis_mm_mu, xlabel, J, M_PT, thinning)
    plot3 = plotParameter(getGamma1(PT_labda, M_PT), "labda", real_param.labda, axis_mm_labda, xlabel, M_PT, thinning)
    plot4J = plotParameters(getGamma1(PT_psi, J, M_PT), "psi", real_param.psi, axis_mm_psi, xlabel, J, M_PT, thinning)
    plot5 = plotParameter(getGamma1(PT_acceptance, M_acceptance), "acceptance", 0, axis_mm_acceptance, xlabel, M_acceptance, thinning, false)
    return fillPlotArray(plot1, plot2J, plot3, plot4J, plot5, J)
end


function plotParametersMH(real_param, MH_labda, MH_psi, MH_mu, MH_tau, MHacceptance, J, M_MH, burnin, axis_mm_labda, axis_mm_psi, axis_mm_mu, axis_mm_tau, axis_mm_acceptance, thinning)
    xlabel = "MH Solo Iteration"
    plot1 = plotParameter(MH_tau, "tau", real_param.tau, axis_mm_tau, xlabel, M_MH, thinning)
    plot2J = plotParameters(reshape(MH_mu, J, M_MH), "mu", real_param.mu, axis_mm_mu, xlabel, J, M_MH, thinning)
    plot3 = plotParameter(MH_labda, "labda", real_param.labda, axis_mm_labda, xlabel, M_MH, thinning)
    plot4J = plotParameters(reshape(MH_psi, J, M_MH), "psi", real_param.psi, axis_mm_psi, xlabel, J, M_MH, thinning)
    plot5 = plotParameter(MHacceptance, "acceptance", 0, axis_mm_acceptance, xlabel, length(MHacceptance), thinning, false)
    return fillPlotArray(plot1, plot2J, plot3, plot4J, plot5, J)
end

function getParam(labda::Array{Float64,1}, psi::Array{Float64,2}, mu::Array{Float64,2}, tau::Array{Float64,1}, J::Int64, M::Int64, burnin::Int64)
    nonBurnIn = (burnin + 1):M
    param_mu = zeros(Float64, J)
    param_psi = zeros(Float64, J)
    for j in 1:J
        param_mu[j] = mean(mu[j, nonBurnIn])
        param_psi[j] = mean(psi[j, nonBurnIn])
    end
    param_labda = mean(labda[nonBurnIn])
    param_tau = mean(tau[nonBurnIn])
    return MHiterationParameters(param_labda, param_psi, param_mu, param_tau, [0.0, 0.0])
end


function plotFDensity(real_param::MHiterationParameters, MH_param::MHiterationParameters, PT_param::MHiterationParameters, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, J, M_MH, M_PT, burnin, burninlabel, colors, thinning)
    xminmax = zeros(2)
    xminmax[1] = min(minimum(real_param.mu) - 3 * sqrt(1 / real_param.tau), minimum(MH_param.mu) - 3 * sqrt(1 / MH_param.tau), minimum(PT_param.mu) - 3 * sqrt(1 / PT_param.tau))
    xminmax[2] = max(maximum(real_param.mu) + 3 * sqrt(1 / real_param.tau), maximum(MH_param.mu) + 3 * sqrt(1 / MH_param.tau), maximum(PT_param.mu) + 3 * sqrt(1 / PT_param.tau))
    xminmax = round(xminmax, 1)

    x = collect(xminmax[1]:1e-2:xminmax[2])
    real_phi = Array(Distribution, 0); real_rho = zeros(Float64, J);
    real_J = length(real_param.mu)
    for j in 1:real_J
        push!(real_phi, Normal(real_param.mu[j], sqrt(1 / real_param.tau)))
        real_rho[j] = real_param.psi[j] / real_param.labda
    end

    y_real = zeros(Float64, length(x)); y_MH = zeros(Float64, length(x)); y_PT = zeros(Float64, length(x))
#     info(size(MH_psi))
#     info(size(PT_psi))
    for i in 1:length(x)
        tmpMH = 0.0;
        tmpPT = 0.0;
        for m in burnin:M_MH
            tmpMH += getDensityf(J, MH_labda[m], MH_psi[:,m], MH_mu[:,m], MH_tau[m], x[i])
        end        
        for m in burnin:M_PT
            tmpPT += getDensityf(J, PT_labda[1, m], PT_psi[:, 1,m], PT_mu[:, 1,m], PT_tau[1, m], x[i])
        end
        y_MH[i] = tmpMH / (M_MH - burnin)
        y_PT[i] = tmpPT / (M_PT - burnin)
        for j in 1:real_J
            y_real[i] += real_rho[j] * pdf(real_phi[j], x[i])
        end
    end
    normdiffMH = norm(y_real - y_MH)
    normdiffPT = norm(y_real - y_PT)
    
    lambdadiffRealMHPT = Array(Float64, 5);
    lambdadiffRealMHPT[1] = real_param.labda
    lambdadiffRealMHPT[2] = MH_param.labda
    lambdadiffRealMHPT[3] = abs(lambdadiffRealMHPT[1] - lambdadiffRealMHPT[2])
    lambdadiffRealMHPT[4] = PT_param.labda
    lambdadiffRealMHPT[5] = abs(lambdadiffRealMHPT[1] - lambdadiffRealMHPT[4])
#     info(lambdadiffRealMHPT)


    myplot = plot(layer(y = y_real, x = x, Geom.line, Theme(line_style = Gadfly.get_stroke_vector(:dot))),
                    layer(y = y_MH, x = x, Geom.line, Theme(default_color = colors[1])),
                    layer(y = y_PT, x = x, Geom.line, Theme(default_color = colors[2])),
                    Guide.xlabel("x"),
                    Guide.ylabel("f(x)"),
                    Coord.Cartesian(xmin = xminmax[1], xmax = xminmax[2]),
                    )
#     append!(myplot.layers, layer(y = y_real, x = x, Geom.line, Theme(line_style = Gadfly.get_stroke_vector(:dot))))
    return myplot, normdiffMH, normdiffPT, lambdadiffRealMHPT
end


function plotFDensities(real_param::MHiterationParameters, PT_param::MHiterationParameters, PT_labda, PT_psi, PT_mu, PT_tau, J, M_PT, burnin, burninlabel, colors, thinning)
    xminmax = zeros(2)
    xminmax[1] = min(minimum(real_param.mu) - 3 * sqrt(1 / real_param.tau), minimum(PT_param.mu) - 3 * sqrt(1 / PT_param.tau))
    xminmax[2] = max(maximum(real_param.mu) + 3 * sqrt(1 / real_param.tau), maximum(PT_param.mu) + 3 * sqrt(1 / PT_param.tau))
    xminmax = round(xminmax, 1)

    x = collect(xminmax[1]:1e-2:xminmax[2])
    real_phi = Array(Distribution, 0); real_rho = zeros(Float64, J);
    real_J = length(real_param.mu)
    for j in 1:real_J
        push!(real_phi, Normal(real_param.mu[j], sqrt(1 / real_param.tau)))
        real_rho[j] = real_param.psi[j] / real_param.labda
    end

    G = length(PT_labda[:, 1])
    y_real = zeros(Float64, length(x)); y_PT = zeros(Float64, G, length(x))

    for i in 1:length(x)
        for g in 1:G
            tmpPT = 0.0;
            for m in burnin:M_PT
                tmpPT += getDensityf(J, PT_labda[g, m], PT_psi[:, g, m], PT_mu[:, g, m], PT_tau[g, m], x[i])
            end
            y_PT[g, i] = tmpPT / (M_PT - burnin)
        end
        for j in 1:real_J
            y_real[i] += real_rho[j] * pdf(real_phi[j], x[i])
        end
    end
    normdiffMH = 1
    normdiffPT = 1
#     norm(y_real - y_PT[1,:])
    
    lambdadiffRealMHPT = Array(Float64, 5);
    lambdadiffRealMHPT[1] = real_param.labda
    lambdadiffRealMHPT[2] = 1
    lambdadiffRealMHPT[3] = abs(lambdadiffRealMHPT[1] - lambdadiffRealMHPT[2])
    lambdadiffRealMHPT[4] = PT_param.labda
    lambdadiffRealMHPT[5] = abs(lambdadiffRealMHPT[1] - lambdadiffRealMHPT[4])
#     info(lambdadiffRealMHPT)


    myplot = plot(layer(y = y_real, x = x, Geom.line, Theme(line_style = Gadfly.get_stroke_vector(:dot))),
#                     layer(y = y_PT, x = x, Geom.line, Theme(default_color = colors[2])),
                    Guide.xlabel("x"),
                    Guide.ylabel("f(x)"),
                    Coord.Cartesian(xmin = xminmax[1], xmax = xminmax[2]),
                    )
#     append!(myplot.layers, layer(y = y_real, x = x, Geom.line, Theme(line_style = Gadfly.get_stroke_vector(:dot))))
    for g in 1:G
        append!(myplot.layers, layer(y = y_PT[g, :], x = x, Geom.line, Theme(default_color = colors[g])))
    end
    return myplot, normdiffMH, normdiffPT, lambdadiffRealMHPT
end



function getDensityf(J, labda, psi, mu, tau, x)
#     phi = Array(Distribution, 0);
#     rho = zeros(Float64, J);
#     push!(phi, Normal(mu[j], sqrt(1 / tau)))
    y = 0
    for j in 1:J
        y += (psi[j] / labda) * pdf(Normal(mu[j], sqrt(1 / tau)), x)
    end
    return y
end

function plotFDensity(real_param::MHiterationParameters, MH_param::MHiterationParameters, PT_param::MHiterationParameters, J::Int64, colors, thinning)
    xminmax = zeros(2)
    xminmax[1] = min(minimum(real_param.mu) - 3 * sqrt(1 / real_param.tau), minimum(MH_param.mu) - 3 * sqrt(1 / MH_param.tau), minimum(PT_param.mu) - 3 * sqrt(1 / PT_param.tau))
    xminmax[2] = max(maximum(real_param.mu) + 3 * sqrt(1 / real_param.tau), maximum(MH_param.mu) + 3 * sqrt(1 / MH_param.tau), maximum(PT_param.mu) + 3 * sqrt(1 / PT_param.tau))
    xminmax = round(xminmax, 1)

    x = collect(xminmax[1]:1e-2:xminmax[2])
    real_phi = Array(Distribution, 0); MH_phi = Array(Distribution, 0); PT_phi = Array(Distribution, 0)
    real_rho = zeros(Float64, J); MH_rho = zeros(Float64, J); PT_rho = zeros(Float64, J)
    for j in 1:J
        push!(MH_phi, Normal(MH_param.mu[j], sqrt(1 / MH_param.tau)))
        push!(PT_phi, Normal(PT_param.mu[j], sqrt(1 / PT_param.tau)))
        MH_rho[j] = MH_param.psi[j] / MH_param.labda
        PT_rho[j] = PT_param.psi[j] / PT_param.labda
    end
    real_J = length(real_param.mu)
    for j in 1:real_J
        push!(real_phi, Normal(real_param.mu[j], sqrt(1 / real_param.tau)))
        real_rho[j] = real_param.psi[j] / real_param.labda
    end

    y_real = zeros(Float64, length(x)); y_MH = zeros(Float64, length(x)); y_PT = zeros(Float64, length(x))
    for i in 1:length(x)
        for j in 1:J
            y_MH[i] += MH_rho[j] * pdf(MH_phi[j], x[i])
            y_PT[i] += PT_rho[j] * pdf(PT_phi[j], x[i])
        end
        for j in 1:real_J
            y_real[i] += real_rho[j] * pdf(real_phi[j], x[i])
        end
    end
    normdiffMH = norm(y_real - y_MH)
    normdiffPT = norm(y_real - y_PT)

    myplot = plot(layer(y = y_real, x = x, Geom.line, Theme(line_style = Gadfly.get_stroke_vector(:dot))),
                    layer(y = y_MH, x = x, Geom.line, Theme(default_color = colors[1])),
                    layer(y = y_PT, x = x, Geom.line, Theme(default_color = colors[2])),
                    Guide.xlabel("x"),
                    Guide.ylabel("f(x)"),
                    Coord.Cartesian(xmin = xminmax[1], xmax = xminmax[2]),
                    )
#     append!(myplot.layers, layer(y = y_real, x = x, Geom.line, Theme(line_style = Gadfly.get_stroke_vector(:dot))))
    return myplot, normdiffMH, normdiffPT
end

function plotF(real_param, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, J, M_MH, M_PT, burnin, burninlabel, colors, thinning)
    MHparam = getParam(MH_labda, MH_psi, MH_mu, MH_tau, J, M_MH, burnin)
    PTparam = getParam(getGamma1(PT_labda, M_PT), getGamma1(PT_psi, J, M_PT), getGamma1(PT_mu, J, M_PT), getGamma1(PT_tau, M_PT), J, M_PT, burnin)    
    
#     METHOD 1
#     USING THE POSTERIOR MEAN
#     plotF, normdiffMH, normdiffPT = plotFDensity(real_param, MHparam, PTparam, J, colors, thinning)
#     METHOD 2
#     USING THE MEAN OF F
    plotF, normdiffMH, normdiffPT, lambdadiffRealMHPT = plotFDensity(real_param, MHparam, PTparam, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, J, M_MH, M_PT, burnin, burninlabel, colors, thinning)
    return plotF, normdiffMH, normdiffPT, lambdadiffRealMHPT
end


function plotF(real_param, PT_labda, PT_psi, PT_mu, PT_tau, J, M_PT, burnin, burninlabel, colors, thinning)
    PTparam = getParam(getGamma1(PT_labda, M_PT), getGamma1(PT_psi, J, M_PT), getGamma1(PT_mu, J, M_PT), getGamma1(PT_tau, M_PT), J, M_PT, burnin)    
    plotF, normdiffMH, normdiffPT, lambdadiffRealMHPT = plotFDensities(real_param, PTparam, PT_labda, PT_psi, PT_mu, PT_tau, J, M_PT, burnin, burninlabel, colors, thinning)
# #     METHOD 1
# #     USING THE POSTERIOR MEAN
# #     plotF, normdiffMH, normdiffPT = plotFDensity(real_param, MHparam, PTparam, J, colors, thinning)
# #     METHOD 2
# #     USING THE MEAN OF F
#     plotF, normdiffMH, normdiffPT, lambdadiffRealMHPT = plotFDensities(real_param, PT_param, PT_labda, PT_psi, PT_mu, PT_tau, J, M_PT, burnin, burninlabel, colors, thinning)
#     return plotF, normdiffMH, normdiffPT, lambdadiffRealMHPT
    return plotF, normdiffMH, normdiffPT, lambdadiffRealMHPT
end

function writePlots(dir, real_param, real_J, M_MH, MMHlabel, M_PT, MPTlabel, J, gamma, UpdateModusMH, noSwap, burnin, burninlabel, InfoMessages, verboseInfoMessages, plotsPTAll, plotsPTGamma1, plotsMH, plotsAC, plotsDensity, plotf, thinning)
    InfoMessages && info("Saving plots to file")
    makeDirectory(dir); dir = joinpath(dir, string("labda=", real_param.labda, ", real-J=", real_J, ", J=", J)); makeDirectory(dir)
    fnameExtension = ".pdf"

    fnameBase = string("Lambda=",round(Int, real_param.labda), "-", UpdateModusMH, "-gamma=", gamma, "-M=", round(Int, MPTlabel / 1000), "k-M_MH=", round(Int, MMHlabel / 1000), "k-ALL")
    noSwap && (fnameBase = string(fnameBase, "_noSwap"))
    savefilePath = fnameNew(fnameBase, fnameExtension, dir)

    fnameBasef = string("Lambda=",round(Int, real_param.labda), "-", UpdateModusMH, "-gamma=", gamma, "-M=", round(Int, MPTlabel / 1000), "k-M_MH=", round(Int, MMHlabel / 1000), "k-burnin=", burninlabel,"-f")
    savefilePathf = fnameNew(fnameBasef, fnameExtension, dir)

    draw(PDF(savefilePathf, 10cm, 10cm), plotf)
    draw(PDF(savefilePath, 50cm, (3 + 2 * J) * 10cm), hstack(vstack(plotsPTAll), vstack(plotsPTGamma1), vstack(plotsMH), vstack(plotsAC), vstack(plotsDensity)))
    InfoMessages && verboseInfoMessages && info("Saving complete")
    filePaths = Array(AbstractString, 0)
    push!(filePaths, savefilePath)
    push!(filePaths, savefilePathf)
    return filePaths
end


function writePlots(dir, real_param, real_J, M_MH, MMHlabel, M_PT, MPTlabel, J, gamma, UpdateModusMH, noSwap, burnin, burninlabel, InfoMessages, verboseInfoMessages, plotsPTAll, plotsDensity, plotf, thinning)
    InfoMessages && info("Saving plots to file")
    makeDirectory(dir); dir = joinpath(dir, string("labda=", real_param.labda, ", real-J=", real_J, ", J=", J)); makeDirectory(dir)
    fnameExtension = ".pdf"

    fnameBase = string("Lambda=",round(Int, real_param.labda), "-", UpdateModusMH, "-gamma=", gamma, "-M=", round(Int, MPTlabel / 1000), "k-M_MH=", round(Int, MMHlabel / 1000), "k-ALL")
    noSwap && (fnameBase = string(fnameBase, "_noSwap"))
    savefilePath = fnameNew(fnameBase, fnameExtension, dir)

    fnameBasef = string("Lambda=",round(Int, real_param.labda), "-", UpdateModusMH, "-gamma=", gamma, "-M=", round(Int, MPTlabel / 1000), "k-M_MH=", round(Int, MMHlabel / 1000), "k-burnin=", burninlabel,"-f")
    savefilePathf = fnameNew(fnameBasef, fnameExtension, dir)

    draw(PDF(savefilePathf, 10cm, 10cm), plotf)
#     draw(PDF(savefilePath, 20cm, (3 + 2 * J) * 10cm), hstack(vstack(plotsPTAll), vstack(plotsDensity)))
    
#     draw(PDF(savefilePath, (3 + 2 * J) * 10cm, 20cm), vstack(hstack(plotsPTAll), hstack(plotsDensity)))
    draw(PDF(savefilePath, (2 + J) * 10cm, 40cm), vstack(vstack(hstack(plotsPTAll[1:(J+2)]), hstack(plotsDensity[1:(J+2)])), vstack(hstack(plotsPTAll[(J+2+1):(3 + 2 * J)]), hstack(plotsDensity[(J+2+1):(3 + 2 * J)]))))
#     draw(PDF(savefilePath, (2 + J) * 10cm, 40cm), title("My awesome data", vstack(vstack(hstack(plotsPTAll[1:(J+2)]), hstack(plotsDensity[1:(J+2)])), vstack(hstack(plotsPTAll[(J+2+1):(3 + 2 * J)]), hstack(plotsDensity[(J+2+1):(3 + 2 * J)])))))
    
    
    
    InfoMessages && verboseInfoMessages && info("Saving complete")
    filePaths = Array(AbstractString, 0)
    push!(filePaths, savefilePath)
#     push!(filePaths, savefilePathf)
    return filePaths
end


function writecsvMH(real_param::MHiterationParameters, real_J::Int64, MHparamArray::Array{ParametersCPPMH.MHiterationParameters,1}, MHacceptance::Array{Float64,1}, PTparamArray::Array{ParametersCPPMH.MHiterationParameters,2}, PTacceptance, N::Int64, J::Int64, I::Int64, G::Int64, M_PT::Int64, M_MH::Int64, gamma::Array{Float64,1}, UpdateModusMH::AbstractString, noSwap::Bool, dir::AbstractString, InfoMessages::Bool, verboseInfoMessages::Bool)
    MH_psi = zeros(J, M_MH); MH_mu = zeros(J, M_MH); MH_tau = zeros(M_MH); MH_labda = zeros(M_MH);
    extractData!(MH_labda, MH_psi, MH_mu, MH_tau, J, M_MH, MHparamArray)
    sortRealParamArray!(real_param)
    mu1 = reshape(MH_mu[1, :], M_MH)
    mu2 = reshape(MH_mu[2, :], M_MH)
#     mu3 = reshape(MH_mu[3, :], M_MH)
#     mu4 = reshape(MH_mu[4, :], M_MH)
    psi1 = reshape(MH_psi[1, :], M_MH)
    psi2 = reshape(MH_psi[2, :], M_MH)
#     psi3 = reshape(MH_psi[3, :], M_MH)
#     psi4 = reshape(MH_psi[4, :], M_MH)

    dir = string("/home/wikash/Downloads/sim-results/J=", J, "/"); makeDirectory(dir)
    fnameResultsBase = string("MH_M=", round(Int, M_MH / 1000), "k_", real_param)
    fnameResultsExtension = ".csv"
    savefilePath = fnameNew(fnameResultsBase, fnameResultsExtension, dir)
    writecsv(savefilePath, zip(mu1, mu2, psi1, psi2, MH_tau, MHacceptance))
#     writecsv(savefilePath, zip(mu1, mu2, mu3, psi1, psi2, psi3, MH_tau, MHacceptance))
#     writecsv(savefilePath, zip(mu1, mu2, mu3, mu4, psi1, psi2, psi3, psi4, MH_tau, MHacceptance))
#     writecsv("/home/wikash/Downloads/csv.txt", [reshape(MH_mu[1, :], M_MH), reshape(MH_mu[2, :], M_MH)])
#     writecsv("/home/wikash/Downloads/csv.txt", [MH_mu[1, :], MH_mu[2, :], MH_mu[3, :], MH_mu[4, :], MH_psi[1, :], MH_psi[2, :], MH_psi[3, :], MH_psi[4, :], MH_tau, MHacceptance])
#     writecsv("/home/wikash/Downloads/csv-acceptance.txt", MHacceptance)
#     writedlm(f, zip(x, y))
#
# psi1 t/m psi4, mu1 t/m mu4 en tau
    return
end

function pdfPlotResults(real_param::MHiterationParameters, real_J::Int64, MHparamArray::Array{ParametersCPPMH.MHiterationParameters,1}, MHacceptance::Array{Float64,1}, PTparamArray::Array{ParametersCPPMH.MHiterationParameters,2}, PTacceptance, N::Int64, J::Int64, I::Int64, G::Int64, MPT::Int64, MMH::Int64, gamma::Array{Float64,1}, burninlabel::Int64, UpdateModusMH::AbstractString, noSwap::Bool, thinning::Bool, dir::AbstractString, InfoMessages::Bool, verboseInfoMessages::Bool)
    InfoMessages && info("Start plotting")
    InfoMessages && (verboseInfoMessages && info("Extracting information"))
    if (thinning)
        M_MH = Int(round(MMH / 10))
        M_PT = Int(round(MPT / 10))
    else
        M_MH = MMH
        M_PT = MPT
    end
    
    PT_psi = zeros(J, G, M_PT); PT_mu = zeros(J, G, M_PT); PT_tau = zeros(G, M_PT); PT_labda = zeros(G, M_PT);
    MH_psi = zeros(J, M_MH); MH_mu = zeros(J, M_MH); MH_tau = zeros(M_MH); MH_labda = zeros(M_MH);
    
    if thinning
        burnin = Int(round(burninlabel / 10))
    else
        burnin = burninlabel
    end
    MMHlabel = MMH
    MPTlabel = MPT

    
    # get data from paramArray, pre-process data
    extractData!(MH_labda, MH_psi, MH_mu, MH_tau, J, M_MH, MHparamArray)
    extractData!(PT_labda, PT_psi, PT_mu, PT_tau, J, G, M_PT, PTparamArray)
    sortRealParamArray!(real_param)
    PTacceptance = PTacceptance[:, PTacceptance[1, :] .!= 0] # remove acceptance = 0 entries (swapping)

    # plotting setup
    const PTlabels = setLabels(gamma, G)
    axis_mm_labda, axis_mm_psi, axis_mm_mu, axis_mm_tau, axis_mm_acceptance = setAxisMinMax(real_param, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau)
    const colorsMHPT = [colorant"black", colorant"red"]
    const colornamesMHPT = ["black", "red"]

    # create plots
    InfoMessages && (verboseInfoMessages && info("Creating plots"))
    plotsPTGamma1 = plotParametersPTGamma1(real_param, PT_labda, PT_psi, PT_mu, PT_tau, PTacceptance, J, M_PT, burnin, axis_mm_labda, axis_mm_psi, axis_mm_mu, axis_mm_tau, axis_mm_acceptance, PTlabels, UpdateModusMH, thinning)
    plotsPTAll = plotParametersPTAllGamma(real_param, PT_labda, PT_psi, PT_mu, PT_tau, PTacceptance, J, G, M_PT, burnin, axis_mm_labda, axis_mm_psi, axis_mm_mu, axis_mm_tau, axis_mm_acceptance, PTlabels, colors, UpdateModusMH, noSwap, thinning)
    
    plotsMH = plotParametersMH(real_param, MH_labda, MH_psi, MH_mu, MH_tau, MHacceptance, J, M_MH, burnin, axis_mm_labda, axis_mm_psi, axis_mm_mu, axis_mm_tau, axis_mm_acceptance, thinning)

    plotsAC, ESSAC = plotAutocorrelationsPTMH(MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, PTacceptance, J, G, M_MH, M_PT, burnin, burninlabel, PTlabels, colorsMHPT, thinning)
    if noSwap
        plotsDensity = plotDensitiesLayered(real_param, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, PTacceptance, J, G, M_MH, M_PT, burnin, burninlabel, PTlabels, colors, colornames, thinning)
    else
        plotsDensity = plotDensitiesMHPT(real_param, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, PTacceptance, J, G, M_MH, M_PT, burnin, burninlabel, PTlabels, colorsMHPT, colornamesMHPT, thinning)
    end
    
    if noSwap
        plotf, normdiffMH, normdiffPT, lambdadiffRealMHPT = plotF(real_param, PT_labda, PT_psi, PT_mu, PT_tau, J, M_PT, burnin, burninlabel, colors, thinning)
    else
        plotf, normdiffMH, normdiffPT, lambdadiffRealMHPT = plotF(real_param, MH_labda, MH_psi, MH_mu, MH_tau, PT_labda, PT_psi, PT_mu, PT_tau, J, M_MH, M_PT, burnin, burninlabel, colorsMHPT, thinning)
    end
    

    # write plots
    if noSwap
        filePaths = writePlots(dir, real_param, real_J, M_MH, MMHlabel, M_PT, MPTlabel, J, gamma, UpdateModusMH, noSwap, burnin, burninlabel, InfoMessages, verboseInfoMessages, plotsPTAll, plotsDensity, plotf, thinning)
    else
        filePaths = writePlots(dir, real_param, real_J, M_MH, MMHlabel, M_PT, MPTlabel, J, gamma, UpdateModusMH, noSwap, burnin, burninlabel, InfoMessages, verboseInfoMessages, plotsPTAll, plotsPTGamma1, plotsMH, plotsAC, plotsDensity, plotf, thinning)
    end
    
    return filePaths, ESSAC, normdiffMH, normdiffPT, lambdadiffRealMHPT
end

end

#, x=collect(10:10:10*length(a))
