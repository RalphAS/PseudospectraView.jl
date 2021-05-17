using PseudospectraView
using QML
using Pseudospectra
const PSA=Pseudospectra
using Test

@testset "Basic" begin
    A = PSA.landau_fox_li(200,16)

    # set the time bomb
    PseudospectraView._doomed[] = true

    psagui(A; backend=:gr)

    app_data = PseudospectraView.PSApp.ps_data
    @test size(app_data.ps_dict[:schur_mtx]) == size(A)
end

