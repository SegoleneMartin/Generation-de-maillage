
workspace()

import Base: convert, ==, isless, +, -, *, sin, cos, angle

typealias RatInt Union{Rational, Integer}  # je crois qu'ici on crée un nouveau type (RatInt) qui est le type "rationnel ou entier"

##########
# ANGLES #
##########

"""
Encodes the angle
    2pi/N*(val + slice)
where val is always assumed to be in [0,1)
"""
immutable Angle{N}              """ Angle est un 'composite type' dont les "coordonnées" sont val et slice. 
    Il prend le paramètre abstrait N qui normalement est un type (par ext N=Int ou N=Float64) (dans ce cas ci a quoi sert-il?)
'immutable' dit qu'on ne peut pas changer les valeurs de val et slice un fois qu'on leur a affecté des valeurs """
    val::Rational{Int}          # value in the camembert slice 
""" The :: operator can be used to give a type to variables in programs. It  declares the variable to always have the specified type
    Ici on dit que val est un rationnel et que slice est un entier"""
    slice::Int                  # camembert slice
end

Angle{T<:RatInt}(val::T, N::Int) = Angle(val, 1, N)
  #Angle{T<:RatInt} declares that Angle takes one type parameter of an RatInt type
#  creates an angle with slice=0 simply omitting the slice parameter (and hence having as parameters only val and N)
"""Recall that in julia the same function has different definitions depending on the number and types of its arguments, this behavior is called "multiple dispatch" """
# N=slice, donc slice = 0
Angle{T<:Integer}(val::T, n::Int, N::Int) = Angle(convert(Rational{Int},val), n,N)
#allows to create angles with integer values (which will hence always result in angles with 0 value and slice+value as slice)


function Angle(val::Rational{Int}, n::Int, N::Int)  """ A constructor is a function which is used to "build" a type, and has usually
the same name as the type itself. Whenever you define a type as above, julia automatically defines a constructor 
    but we want something that can enforce the above bounds val in [0,1[ and slice<N. On met N en argument pour que ce soit plus pratique"""
    if 0<=val<1
        Angle{N}(val, n)
    elseif val>=0
        Angle{N}(val%1, (div(val,1)+n)%N)
    else                        # val<0 
        Angle{N}(1+val%1, (N+div(val,1)+n-1)%N)
    end
end  #Essentially, we check if val is in  [0,1) and if not we shift its integer part to slice. Ainsi, contrairement au type "Angle", le constructor
# Angle permet de prendre l'angle modulo 2pi

camembert{N}(::Angle{N}) = N  #Angle{N} est un composite type, camembert prendre des arguments de type Angle{N}
slice(x::Angle) = x.slice   #slice prend Angle en argument et rend la coordonée "slice" de Angle
rotate(a::Angle, n::Int) = Angle{camembert(a)}(a.val, (a.slice + n)%camembert(a)) # rotation prend un angle en argument et renvoie l'angle + n parts
# du camembert

convert(T::Type{AbstractFloat}, x::Angle) = (convert(T, x.val)+x.slice)*2pi/camembert(x) # ??
#Cette première ligne est-elle vraiment utile?
convert{T<:AbstractFloat}(::Type{T}, x::Angle) = (convert(T, x.val)+x.slice)*2pi/camembert(x) # pour tous les sous types de AbstracFloat
convert{T<:Rational{Int}}(::Type{Angle}, x::T) = Angle(x, 1)   #Que font ces deux lignes?
convert{T<:Integer}(::Type{Angle}, x::T) = Angle(x//1, 1)

==(a::Angle, b::Angle) = camembert(a) == camembert(b) && (a.val == b.val && a.slice == b.slice) # on définit l'égalité de deux angles a et b
#PB: A aucun moment dans cette definition on fait le modulo. --> Angle{6}(1,3)==Angle{6}(0,4) renvoie false alors que ce sont les mêmes angles
isless(a::Angle, b::Angle) = (a.val + a.slice) < (b.val + b.slice)  #on définit ce qu'est un plus petit angle

function +(a::Angle, b::Angle)   # on crée une fonction qui additionne deux angles
    if camembert(a) == camembert(b)
        Angle(a.val+b.val, a.slice + b.slice, camembert(a))
    else
        error("Angles needs to have the same camembert slices")
    end
end

+{T<:RatInt}(a::Angle, b::T) = Angle(a.val+b, a.slice, camembert(a)) # ATTENTION! 
+{T<:RatInt}(b::T, a::Angle) = a+b

-(a::Angle) = Angle(-a.val, -a.slice, camembert(a))  # définit l'opposé d'un angle
-(a::Angle, b::Angle) = a+(-b)    # définit la soustraction
-{T<:RatInt}(a::Angle, b::T) = -a+b
-{T<:RatInt}(b::T, a::Angle) = b+(-a)

*{T<:RatInt}(a::Angle, b::T) = Angle(a.val*b, a.slice*b, camembert(a))
*{T<:RatInt}(b::T, a::Angle) = a*b

sin(a::Angle) = sin(convert(Float64, a))
cos(a::Angle) = cos(convert(Float64, a))


###############
# FREQUENCIES #
###############

"""
2D frequency in polar coordinates
"""
immutable Frequency{N, T<:Real}      
    λ::T      # module de l'angle
    ω::Angle{N}
end

Frequency{R<:Real}(a::R, x...) = Frequency{x[end], R}(a, Angle(x...))   # il faut voir x... comme un tuple (val, slice, N). x[end]=N
# en gros Frequency(2.5, Angle{6}(1//2,1)) rend la même chose que Frequency(2.5, 1//2, 1, 6)

slice(x::Frequency) = slice(x.ω)    # on définit slice lorsque slice prend en argument une fréquence
camembert{T<:Real,N}(::Frequency{N,T}) = N
Base.angle(x::Frequency) = convert(Float64, x.ω)
rotate(x::Frequency, n::Int) = Frequency(x.λ, rotate(x.ω,n))   # cf rotate précedent

==(x::Frequency, y::Frequency) = x.λ == y.λ && x.ω == y.ω
isless(x::Frequency, y::Frequency) = x.λ < y.λ || (x.λ == y.λ && x.ω < y.ω)

-(a::Frequency) = Frequency(a.λ, -a.ω)

"""
Product of two frequencies.
"""
composition{T<:Real,N}(x::Frequency{N,T}, y::Frequency{N,T}) = Frequency{N,T}(x.λ*y.λ, x.ω-y.ω)  

"""
Cartesian coordinates of the frequency.
"""
cart(x::Frequency) = (x.λ*cos(x.ω), x.λ*sin(x.ω))
cart{T<:Real, N}(v::Vector{Frequency{N, T}}) = (T[x.λ*cos(angle(x)) for x in v], T[x.λ*sin(angle(x)) for x in v])  # array 2xN. Vector(A,B, ...N) avec A, B...N des fréquences

normalize(x::Frequency) = rotate(x, -slice(x))  # ramène à la slice 0

###################
# BISPECTRAL SETS #
###################

"""
A bispectral set is a special vector of frequencies.
"""
immutable BispectralSet{N, T<:Real} <: AbstractArray{Frequency{N,T},1}  
    pts::Vector{Frequency{N, T}}
end

function BispectralSet{T<:Real, N}(pts::Vector{Frequency{N, T}})   # pourquoi pts n'est pas un vecteur?
    @assert all(x->x<=1,[slice(x) for x in pts]) "Frequencies must be in [0,2π/N)"      # en fait je crois que c'est correct: on veut une fréquence
                                      # entre 0 et 2pi/N, ie un angle dont la slice est 0. Par contre faudrait peut-être remplacer x->x<=1 par x->x=1
#  @assert all(x->x==N,[camembert(x) for x in pts]) "Number of slices must be the N for all frequencies"
    x = issorted(pts) ? pts : sort(pts)  # Tests if a vector is in sorted order and sort it if not (cf def isless pour les fréquences)
    BispectralSet{camembert(pts[1]), T}(x)    # on prendre camembert(pts[1]) puisque de toute façon toutes les freq contenues dans pts ont le même camembert
    end  # en gros cette fonction dit que quand on veut le Bispectralset d'un vector pts, on nous renvoie le BispectralSet du vecteur pts trié

# BispectralSet{T<:Real, N}(N::Int,pts::Vector{Frequency{T}}) = BispectralSet{T}(N,pts)  #???

camembert{N, T<:Real}(::BispectralSet{N,T}) = N

Base.size(E::BispectralSet) = (length(E.pts), camembert(E))   # définition dindon ddDDdd
Base.size(E::BispectralSet, n) = size(E)[n]   # n est l'index (1 ou 2 ?)
# Base.linearindexing(::Type{BispectralSet}) = Base.LinearSlow()
function Base.getindex(E::BispectralSet, i::Int) 
    if 1<= i <= size(E,1) 
        E.pts[i] 
    else
        r = mod1(i, size(E,1))
        d = div(i-r, size(E,1))  # cimer a quoi ca sert le '-r'?
        E[r, d+1]  # index de la matrice
    end
end
Base.getindex(E::BispectralSet, i::Int, n::Int) = 1<=n<=camembert(E) ? rotate(E.pts[i], n-1) : BoundsError(E, [i,n])
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, ::Colon, ns) = Frequency{N,T}[E[i,n] for i in 1:size(E,1), n in ns]   #Colon = tous les vecteurs
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, is, ns) = Frequency{N,T}[E[i,n] for i in is, n in ns]
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, ::Colon) = vec(E[1:size(E,1), 1:size(E,2)])

Base.start(::BispectralSet) = 1
Base.next(E::BispectralSet, state) = (E[state], state+1)
Base.done(E::BispectralSet, s) = s > prod(size(E))   # nombre d'elt

cart(E::BispectralSet) = cart(E[:])  # met en cartesienne les frequences de E

"""
Generates a BispectralSet given a vector (`cutoff`) 
of tuples of the type (n, range), where range are the radii of the frequencies
and n the number of equispaced angles in the camembert slice for each one of them.
"""
function BispectralSet{T<:Range}(N::Int, cutoff::Vector{Tuple{Int64,T}})
    fill_freqs(x) = fill_freq_vector(N, x...)
    v = mapreduce(fill_freqs, append!, cutoff)   # pour tout elt dans cutoff on fait fill_freq --> vecteur de vecteurs de freq puis on fait reduce avec append (colle)
    BispectralSet(v)
end

"""
Auxiliary function to fill the frequency vector.
"""
function fill_freq_vector(N::Int, angles::Int, radii)  # cutoff=(angles, radii)
    n = length(radii)
    v = Vector{Frequency{N, Float64}}(n*angles)
    pos = 1
    for i in 1:n, j = 0:angles-1
        @inbounds v[pos] = Frequency(radii[i], j//angles, N) #  j//angles est la val, N est le camember. On omet slice --> slice=1   # on ne vérifie pas que i est bien compris entre 1 et n
        pos += 1
    end
    v
end

"""
Simple search by bisection
"""
function Base.findin(x::Frequency, E::BispectralSet)  # cherche si une freq est dans l'ensemble Bispectral et si oui renvoie l'index
    y = normalize(x)
    left = 1
    (y == E[left]) && (return left, slice(x)+1)
    right = size(E,1)
    (y == E[right]) && (return right, slice(x)+1)
    middle = round(Int, size(E,1)/2)
    for i in 1:size(E,1)
        ((middle <= left) || middle >= right) && break
        if y == E[middle]
            return middle, slice(x)+1 
        elseif    y > E[middle]
            left = middle
        else                    #if y<E[middle]
            right = middle
        end
        middle = round(Int, (right+left)/2)
    end
    return 0,0
end

"""
Extracts a vector `angles` such that `angles[i]` contains the
number of elements in the camembert slice at the radius `i`.
"""
function camembert_angles{N,T<:Real}(E::BispectralSet{N,T})
    ρs = radii(E)
    cam_angles = angles_def(E)
    angles = zeros(Int, length(ρs))
    pos_ρ = 1
    for i in 1:length(ρs)
        while ρs[i] == E[pos_ρ].λ
            angles[i] += 1
            pos_ρ+=1
        end
    end
    angles
end

#############################
# EXTRACT RADIAL DEFINITION #
#############################

Base.den{T<:Integer}(x::Vector{Rational{T}}) = map(den, x)
Base.angle(E::BispectralSet) = Rational{Int}[E[i].ω.val for i in 1:size(E,1)] |> unique
radii{T<:Real,N}(E::BispectralSet{N,T}) = T[E[i].λ for i in 1:size(E,1)] |> unique
angles_def(E::BispectralSet) = angle(E) |> den |> lcm 

function Base.gcd{T<:Integer}(a::Rational{T}, b::Rational{T})
    l = lcm(den(a),den(b))
    gcd(num(a)*gcd(den(a),l), num(b)*gcd(den(b),l))//l
end
Base.gcd{T<:Real}(a::T, b::T) = gcd(rationalize(Int,a), rationalize(Int, b))
Base.gcd{T<:Real, R<:Integer}(a::Rational{R}, b::T) = gcd(a, rationalize(Int, b))

radius_def(E::BispectralSet) =  radius_def(radii(E))
function radius_def{T<:Real}(x::Vector{T})   # trouve le plus petit ecart dans radii
    v = (x-circshift(x,1))[2:end]  # copie du vecteur décalé de 1
    foldl(gcd, v)   # reduce from the left
end


