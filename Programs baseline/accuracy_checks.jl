# Loop through the keys in dict1 and compare with dict2
function percentage_difference(value1, value2)
    return 100 * (value2 - value1) / value1
end


function dict_comparison(ss::Dict, ss2::Dict)
    for key in keys(ss)
        if haskey(ss2, key)
            value1 = ss[key]
            value2 = ss2[key]
            percent_diff = percentage_difference(value1, value2)
            if percent_diff > 1e-4
                println("Key: $key, Dict1: $value1, Dict2: $value2, Percentage Difference: $percent_diff%")
            end
        else
            println("Key $key not found in dict2.")
        end
    end
end 