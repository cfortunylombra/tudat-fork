/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "acceleration.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Convert `AvailableAcceleration` to `json`.
void to_json( json& jsonObject, const AvailableAcceleration& availableAcceleration )
{
    jsonObject = json_interface::stringFromEnum( availableAcceleration, availableAccelerations );
}

//! Convert `json` to `AvailableAcceleration`.
void from_json( const json& jsonObject, AvailableAcceleration& availableAcceleration )
{
    availableAcceleration = json_interface::enumFromString( jsonObject.get< std::string >( ), availableAccelerations );
}

} // namespace basic_astrodynamics


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AccelerationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< AccelerationSettings >& accelerationSettings )
{
    if ( accelerationSettings )
    {
        using namespace basic_astrodynamics;
        using namespace json_interface;
        using K = Keys::Acceleration;

        const AvailableAcceleration accelerationType = accelerationSettings->accelerationType_;

        switch ( accelerationType ) {
        case undefined_acceleration:
        case aerodynamic:
        case cannon_ball_radiation_pressure:
        {
            jsonObject[ K::type ] = accelerationType;
            return;
        }
        case point_mass_gravity:
        case third_body_point_mass_gravity:
        {
            jsonObject[ K::type ] = point_mass_gravity;
            return;
        }
        case spherical_harmonic_gravity:
        case third_body_spherical_harmonic_gravity:
        {
            jsonObject[ K::type ] = spherical_harmonic_gravity;
            boost::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicAccelerationSettings
                    = boost::dynamic_pointer_cast< SphericalHarmonicAccelerationSettings >( accelerationSettings );
            jsonObject[ K::maximumDegree ] = sphericalHarmonicAccelerationSettings->maximumDegree_;
            jsonObject[ K::maximumOrder ] = sphericalHarmonicAccelerationSettings->maximumOrder_;
            return;
        }
        case mutual_spherical_harmonic_gravity:
        case third_body_mutual_spherical_harmonic_gravity:
        {
            jsonObject[ K::type ] = mutual_spherical_harmonic_gravity;
            boost::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicAccelerationSettings
                    = boost::dynamic_pointer_cast< MutualSphericalHarmonicAccelerationSettings >(
                        accelerationSettings );
            jsonObject[ K::maximumDegreeOfBodyExertingAcceleration ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfBodyExertingAcceleration_;
            jsonObject[ K::maximumOrderOfBodyExertingAcceleration ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumOrderOfBodyExertingAcceleration_;
            jsonObject[ K::maximumDegreeOfBodyUndergoingAcceleration ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfBodyUndergoingAcceleration_;
            jsonObject[ K::maximumOrderOfBodyUndergoingAcceleration ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumOrderOfBodyUndergoingAcceleration_;
            jsonObject[ K::maximumDegreeOfCentralBody ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumDegreeOfCentralBody_;
            jsonObject[ K::maximumOrderOfCentralBody ] =
                    mutualSphericalHarmonicAccelerationSettings->maximumOrderOfCentralBody_;
            return;
        }
        case relativistic_correction_acceleration:
        {
            boost::shared_ptr< RelativisticAccelerationCorrectionSettings > relativisticAccelerationCorrectionSettings
                    = boost::dynamic_pointer_cast< RelativisticAccelerationCorrectionSettings >( accelerationSettings );
            jsonObject[ K::calculateSchwarzschildCorrection ] =
                    relativisticAccelerationCorrectionSettings->calculateSchwarzschildCorrection_;
            jsonObject[ K::calculateLenseThirringCorrection ] =
                    relativisticAccelerationCorrectionSettings->calculateLenseThirringCorrection_;
            jsonObject[ K::calculateDeSitterCorrection ] =
                    relativisticAccelerationCorrectionSettings->calculateDeSitterCorrection_;
            jsonObject[ K::primaryBody ] =
                    relativisticAccelerationCorrectionSettings->primaryBody_;
            jsonObject[ K::centralBodyAngularMomentum ] =
                    relativisticAccelerationCorrectionSettings->centralBodyAngularMomentum_;
            return;
        }
        case empirical_acceleration:
        {
            boost::shared_ptr< EmpiricalAccelerationSettings > empiricalAccelerationSettings
                    = boost::dynamic_pointer_cast< EmpiricalAccelerationSettings >( accelerationSettings );
            jsonObject[ K::constantAcceleration ] = empiricalAccelerationSettings->constantAcceleration_;
            jsonObject[ K::sineAcceleration ] = empiricalAccelerationSettings->sineAcceleration_;
            jsonObject[ K::cosineAcceleration ] = empiricalAccelerationSettings->cosineAcceleration_;
            return;
        }
        default:
            throw std::runtime_error( stringFromEnum( accelerationType, availableAccelerations )
                                      + " not supported by json_interface." );
        }
    }
}

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< AccelerationSettings >& accelerationSettings )
{
    using namespace json_interface;
    using namespace basic_astrodynamics;
    using K = Keys::Acceleration;

    // Get acceleration type
    const AvailableAcceleration accelerationType = getValue< AvailableAcceleration >( jsonObject, K::type );

    if ( accelerationType == third_body_point_mass_gravity ||
         accelerationType == third_body_spherical_harmonic_gravity ||
         accelerationType == third_body_mutual_spherical_harmonic_gravity )
    {
        std::cerr << "Whether a body will cause a third-body acceleration is determined internally "
                  << "by Tudat based on the propagation settings.\nRemove \"thirdBody\" from \""
                  << stringFromEnum( accelerationType, availableAccelerations )
                  << "\" at key " << getKeyPath( jsonObject )
                  << " to silence this warning." << std::endl;
    }

    switch ( accelerationType ) {
    case undefined_acceleration:
    case aerodynamic:
    case cannon_ball_radiation_pressure:
    {
        accelerationSettings = boost::make_shared< AccelerationSettings >( accelerationType );
        return;
    }
    case point_mass_gravity:
    case third_body_point_mass_gravity:
    {
        accelerationSettings = boost::make_shared< AccelerationSettings >( point_mass_gravity );
        return;
    }
    case spherical_harmonic_gravity:
    case third_body_spherical_harmonic_gravity:
    {
        accelerationSettings = boost::make_shared< SphericalHarmonicAccelerationSettings >(
                    getNumeric< int >( jsonObject, K::maximumDegree ),
                    getNumeric< int >( jsonObject, K::maximumOrder ) );
        return;
    }
    case mutual_spherical_harmonic_gravity:
    case third_body_mutual_spherical_harmonic_gravity:
    {
        MutualSphericalHarmonicAccelerationSettings defaults( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
        accelerationSettings = boost::make_shared< MutualSphericalHarmonicAccelerationSettings >(
                    getNumeric< int >( jsonObject, K::maximumDegreeOfBodyExertingAcceleration ),
                    getNumeric< int >( jsonObject, K::maximumOrderOfBodyExertingAcceleration ),
                    getNumeric< int >( jsonObject, K::maximumDegreeOfBodyUndergoingAcceleration ),
                    getNumeric< int >( jsonObject, K::maximumOrderOfBodyUndergoingAcceleration ),
                    getNumeric( jsonObject, K::maximumDegreeOfCentralBody, defaults.maximumDegreeOfCentralBody_ ),
                    getNumeric( jsonObject, K::maximumOrderOfCentralBody, defaults.maximumOrderOfCentralBody_ ) );
        return;
    }
    case relativistic_correction_acceleration:
    {
        RelativisticAccelerationCorrectionSettings defaults;
        accelerationSettings = boost::make_shared< RelativisticAccelerationCorrectionSettings >(
                    getValue( jsonObject, K::calculateSchwarzschildCorrection,
                              defaults.calculateSchwarzschildCorrection_ ),
                    getValue( jsonObject, K::calculateLenseThirringCorrection,
                              defaults.calculateLenseThirringCorrection_ ),
                    getValue( jsonObject, K::calculateDeSitterCorrection, defaults.calculateDeSitterCorrection_ ),
                    getValue( jsonObject, K::primaryBody, defaults.primaryBody_ ),
                    getValue( jsonObject, K::centralBodyAngularMomentum, defaults.centralBodyAngularMomentum_ ) );
        return;
    }
    case empirical_acceleration:
    {
        EmpiricalAccelerationSettings defaults;
        accelerationSettings = boost::make_shared< EmpiricalAccelerationSettings >(
                    getValue( jsonObject, K::constantAcceleration, defaults.constantAcceleration_ ),
                    getValue( jsonObject, K::sineAcceleration, defaults.sineAcceleration_ ),
                    getValue( jsonObject, K::cosineAcceleration, defaults.cosineAcceleration_ ) );
        return;
    }
    default:
        throw std::runtime_error( stringFromEnum( accelerationType, availableAccelerations )
                                  + " not supported by json_interface." );
    }
}

} // namespace simulation_setup

} // namespace tudat
