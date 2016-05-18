# gobblegobblee

import Base: convert, ==, isless, +, -, *, sin, cos, angle

typealias RatInt Union{Rational, Integer}

##########
# ANGLES #
##########

"""
Encodes the angle
    2pi/N*(val + slice)
where val is always assumed to be in [0,1)
"""
immutable Angle{N}              # declares that Angle parameters are always of the same type
    val::Rational{Int}          # value in the camembert slice # The :: operator can be used to attach type annotations to expressions 
                                # and variables in programs. Ici on dit que val est un rationnel
    slice::Int                  # camembert slice number
end

Angle{T<:RatInt}(val::T, N::Int) = Angle(val, 1, N) #pour remettre à 0 une fois le tour du camembert fait #ATTENTION 
Angle{T<:Integer}(val::T, n::Int, N::Int) = Angle(convert(Rational{Int},val), n, N) # trois paramètres wtf
function Angle(val::Rational{Int}, n::Int, N::Int) #pourquoi là trois paramètres et ligne 39 y en a que 2?
    if 0<=val<1
        Angle{N}(val, n)
    elseif val>=0
        Angle{N}(val%1, (div(val,1)+n)%N) #pourquoi on ajoute div?
    else                        # val<0 
        Angle{N}(1+val%1, (N+div(val,1)+n-1)%N) #pourquoi y a un N?
    end
end

camembert{N}(::Angle{N}) = N #camembert prend des éléments de type angle (N) + nbr de slices
slice(x::Angle) = x.slice #slice prend un angle en paramètre et ressort x.slice + pas la mm chose que celui d'avant
rotate(a::Angle, n::Int) = Angle{camembert(a)}(a.val, (a.slice + n)%camembert(a)) #rotation de n slice

convert(T::Type{AbstractFloat}, x::Angle) = (convert(T, x.val)+x.slice)*2pi/camembert(x) #met l'angle en radians
convert{T<:AbstractFloat}(::Type{T}, x::Angle) = (convert(T, x.val)+x.slice)*2pi/camembert(x) #pas très utile
convert{T<:Rational{Int}}(::Type{Angle}, x::T) = Angle(x, 1) #angle de 0 à 2pi
convert{T<:Integer}(::Type{Angle}, x::T) = Angle(x//1, 1) #same, + transforme en rationnel

==(a::Angle, b::Angle) = camembert(a) == camembert(b) && (a.val == b.val && a.slice == b.slice) #== fonction, on définit l'égalité, && et, || ou
isless(a::Angle, b::Angle) = (a.val + a.slice) < (b.val + b.slice) #inégalité entre les angles

function +(a::Angle, b::Angle)
    if camembert(a) == camembert(b)
        Angle(a.val+b.val, a.slice + b.slice, camembert(a))
    else
        error("Angles needs to have the same camembert slices")
    end
end

+{T<:RatInt}(a::Angle, b::T) = Angle(a.val+b, a.slice, camembert(a)) # ATTENTION! 
+{T<:RatInt}(b::T, a::Angle) = a+b 

-(a::Angle) = Angle(-a.val, -a.slice, camembert(a)) #opposé, pas soustraction
-(a::Angle, b::Angle) = a+(-b)
-{T<:RatInt}(a::Angle, b::T) = -a+b
-{T<:RatInt}(b::T, a::Angle) = b+(-a)

*{T<:RatInt}(a::Angle, b::T) = Angle(a.val*b, a.slice*b, camembert(a))
*{T<:RatInt}(b::T, a::Angle) = a*b

sin(a::Angle) = sin(convert(Float64, a))#donc convert définit bien en radian
cos(a::Angle) = cos(convert(Float64, a))


###############
# FREQUENCIES #
###############

"""
2D frequency in polar coordinates
"""
immutable Frequency{N, T<:Real} #lamda=module
    λ::T
    ω::Angle{N}
end

Frequency{R<:Real}(a::R, x...) = Frequency{x[end], R}(a, Angle(x...)) #what type is this?

slice(x::Frequency) = slice(x.ω)
camembert{T<:Real,N}(::Frequency{N,T}) = N
Base.angle(x::Frequency) = convert(Float64, x.ω) #Base:module pcq angle - définit avant
rotate(x::Frequency, n::Int) = Frequency(x.λ, rotate(x.ω,n)) #ressort la fréquence après rotation - comme avant

==(x::Frequency, y::Frequency) = x.λ == y.λ && x.ω == y.ω 
isless(x::Frequency, y::Frequency) = x.λ < y.λ || (x.λ == y.λ && x.ω < y.ω)

-(a::Frequency) = Frequency(a.λ, -a.ω)

"""
Product of two frequencies.
"""
composition{T<:Real,N}(x::Frequency{N,T}, y::Frequency{N,T}) = Frequency{N,T}(x.λ*y.λ, x.ω-y.ω) #exp(produit scalaire)

"""
Cartesian coordinates of the frequency.
"""
cart(x::Frequency) = (x.λ*cos(x.ω), x.λ*sin(x.ω))
    cart{T<:Real, N}(v::Vector{Frequency{N, T}}) = (T[x.λ*cos(angle(x)) for x in v], T[x.λ*sin(angle(x)) for x in v]) #vect : une liste. on créé une matrice en rentrant plusieurs valeurs - array 2xn coordonées cart

normalize(x::Frequency) = rotate(x, -slice(x)) #ramène à la slice 0

###################
# BISPECTRAL SETS #
###################

"""
A bispectral set is a special vector of frequencies.
"""
immutable BispectralSet{N, T<:Real} <: AbstractArray{Frequency{N,T},1}
    pts::Vector{Frequency{N, T}}
end

function BispectralSet{T<:Real, N}(pts::Vector{Frequency{N, T}}) 
    @assert all(x->x<=1,[slice(x) for x in pts]) "Frequencies must be in [0,2π/N)" #vérifie que slice=0 --> fréquence dans 2pi/N
#     @assert all(x->x==N,[camembert(x) for x in pts]) "Number of slices must be the N for all frequencies"
    x = issorted(pts) ? pts : sort(pts) # garanti que c'est toujours trié
    BispectralSet{camembert(pts[1]), T}(x)
end

# BispectralSet{T<:Real, N}(N::Int,pts::Vector{Frequency{T}}) = BispectralSet{T}(N,pts)

camembert{N, T<:Real}(::BispectralSet{N,T}) = N

Base.size(E::BispectralSet) = (length(E.pts), camembert(E))
Base.size(E::BispectralSet, n) = size(E)[n]
# Base.linearindexing(::Type{BispectralSet}) = Base.LinearSlow()
function Base.getindex(E::BispectralSet, i::Int) 
    if 1<= i <= size(E,1) 
        E.pts[i] 
    else
        r = mod1(i, size(E,1))
        d = div(i-r, size(E,1))
        E[r, d+1]
    end
end
Base.getindex(E::BispectralSet, i::Int, n::Int) = 1<=n<=camembert(E) ? rotate(E.pts[i], n-1) : BoundsError(E, [i,n])# E[i,n], si n trop grand: erreur
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, ::Colon, ns) = Frequency{N,T}[E[i,n] for i in 1:size(E,1), n in ns] # ns: plusieurs n, + colon : permet de faire toute la range - one st obligé de def 
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, is, ns) = Frequency{N,T}[E[i,n] for i in is, n in ns] # plusieurs i plusieurs n
Base.getindex{T<:Real,N}(E::BispectralSet{N,T}, ::Colon) = vec(E[1:size(E,1), 1:size(E,2)]) # E[:] - affiche toutes les fréquenes (toutes les rotations sous formes de vecteur)

Base.start(::BispectralSet) = 1 # pour pouvoir créer des itérations. start: premier index
Base.next(E::BispectralSet, state) = (E[state], state+1) # index du prochain élément
Base.done(E::BispectralSet, s) = s > prod(size(E)) # quand on finit ; qd plus grand que le produit de la dim (nombre de fréquences + nbrs de slice - on multplie les deux)

cart(E::BispectralSet) = cart(E[:]) # met en cart toutes les fréquences de E

"""
Generates a BispectralSet given a vector (`cutoff`) 
of tuples of the type (n, range), where range are the radii of the frequencies #radii: liste d'entiers
and n the number of equispaced angles in the camembert slice for each one of them.
"""
function BispectralSet{T<:Range}(N::Int, cutoff::Vector{Tuple{Int64,T}})
    fill_freqs(x) = fill_freq_vector(N, x...) # si x tuple, décompose dans les éléments
    v = mapreduce(fill_freqs, append!, cutoff) # on map (effectue) fill_freqs pour tout elt dans cutoff --> vecteur de vecteurs de fréquences, puis on enchaine tout les vecteurs qu'on a trouvé - reduce avec append!
    BispectralSet(v)
end

"""
Auxiliary function to fill the frequency vector.
""" 
function fill_freq_vector(N::Int, angles::Int, radii) # cutoff=(angles, radii)
    n = length(radii)
    v = Vector{Frequency{N, Float64}}(n*angles)
    pos = 1
    for i in 1:n, j = 0:angles-1
        @inbounds v[pos] = Frequency(radii[i], j//angles, N) # j*2pi/N/angles = value + inbounds - sûr qu'on va pas sortir du vect - pas besoin de vérifier 
        pos += 1
    end
    v
end

"""
Simple search by bisection
"""
function Base.findin(x::Frequency, E::BispectralSet) # donne l'index de la fréquence, recheche par bissection
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
        while ρs[i] == E[pos_ρ].λ #nbr d'éléments dans la slice au rayon i
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
function radius_def{T<:Real}(x::Vector{T}) # on cherche le plus petit écart dans radii
    v = (x-circshift(x,1))[2:end] #copie du vecteur décalé de 1
    foldl(gcd, v) #fold from the left. comme reduce mais commence par la gauche (ordre)
end


